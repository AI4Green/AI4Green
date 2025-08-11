import hashlib
import io
import zipfile
from datetime import datetime, timedelta
from typing import Callable

import pandas as pd
import pytz
from azure.storage.blob import (
    BlobClient,
    BlobSasPermissions,
    BlobServiceClient,
    generate_blob_sas,
)
from flask import abort, current_app
from flask.ctx import AppContext
from sources import models, services
from sources.extensions import db


def make_sas_link(data_export_request: models.DataExportRequest) -> str:
    """
    Generates a Shared Access Signature (SAS) token for accessing a blob in an Azure Blob Storage container,
    and constructs the SAS URL for accessing the blob with the given permissions and validity period.

    Args:
        data_export_request (models.DataExportRequest): An instance of DataExportRequest model relating to the
        data export blob

    Returns:
        str: The SAS URL for accessing the blob with the generated SAS token appended as query parameters.
    """
    blob_service_client = services.file_attachments.connect_to_azure_blob_service()
    blob_name = data_export_request.uuid + get_export_file_extension(
        data_export_request.data_format.value
    )
    blob_client = blob_service_client.get_blob_client(
        container="export-outputs",
        blob=blob_name,
    )
    # Set the start time and expiry time for the SAS token
    start_time = datetime.utcnow()
    expiry_time = datetime.utcnow() + timedelta(minutes=2)

    # Define permissions for the SAS token
    permissions = BlobSasPermissions(read=True)
    # ELN file name should match the root folder name not the export uuid
    if data_export_request.data_format.value == "ELN":
        filename = (
            services.data_export.eln_file.get_root_name(data_export_request) + ".eln"
        )
    else:
        filename = blob_name

    # Generate the SAS token
    sas_token = generate_blob_sas(
        blob_client.account_name,
        blob_client.container_name,
        blob_client.blob_name,
        account_key=blob_service_client.credential.account_key,
        permission=permissions,
        expiry=expiry_time,
        start=start_time,
        content_disposition=f'attachment; filename="{filename}"',
    )

    # Construct the SAS URL
    return f"{blob_client.url}?{sas_token}"


def get_export_file_extension(
    export_format: models.data_export_request.ExportFormat,
) -> str:
    """
    Uses the dictionary to return the file extension which matches the export format, including .zip.

    Args:
        export_format: String enum value from an entry in data_export_request Table

    Returns:
        file extension as a string e.g., .csv
    """

    return {
        "ELN": ".eln",
        "RDF": ".zip",
        "RXN": ".zip",
        "JSON": ".zip",
        "CSV": ".csv",
        "SURF": ".txt",
        "PDF": ".zip",
    }[export_format]


def initiate(current_app_context: AppContext, data_export_request_id: int):
    """
    Starts the process of creating the data export once it has been approved.
    This is called from here to enable working within a threaded process.

    Args:
        current_app_context - the app context required to access the database within a threaded process
        data_export_request_id - the primary key id of the data export request we are making the export file for
    """
    with current_app_context:
        export = DataExport(data_export_request_id)
        export.create()


class DataExport:
    """
    Class to control the creation the data export, once approved.
    If we make a zip file, the intermediate container where we save the files used to make the zip is deleted after use.
    Files in the export-container will be deleted after 7 days. Controlled by Azure Blob Storage Lifecycle policies.
    """

    def __init__(self, data_export_request_id: int):
        """
        Create an instance of the DataExport class and setup class variables.

        Args:
            data_export_request_id - the primary key id of the data export request being exported.
        """
        self.data_export_request = models.DataExportRequest.query.get(
            data_export_request_id
        )
        self.container_name = self.data_export_request.uuid
        self.blob_service_client = (
            services.file_attachments.connect_to_azure_blob_service()
        )

    def create(self):
        """
        Main function called for the class where the appropriate export function is identified and called.
        Once the data export has been made, the requestor is notified.
        """
        export_function = self._get_export_function()
        export_function()
        self._update_requestor()

    def _make_zip(self):
        """Make the zip file and generate an md5 checksum to check file integrity"""
        zip_stream = self._make_zip_file()
        checksum = hashlib.md5(zip_stream.getvalue()).hexdigest()
        # upload file and confirm the upload was successful
        blob_client = self._upload_zip(zip_stream)
        self._confirm_upload(blob_client, checksum)
        # update database entry with hash
        self._update_db_with_hash(checksum)
        # delete container
        self.blob_service_client.delete_container(self.container_name)

    def _get_export_function(self) -> Callable:
        """Returns the function which matches the requested data export format"""
        export_function_dict = {
            "RXN": self._make_rxn_export,
            "RDF": self._make_rdf_export,
            "CSV": self._make_csv_export,
            "JSON": self._make_json_export,
            "SURF": self._make_surf_export,
            "PDF": self._make_pdf_export,
            "ELN": self._make_eln_export,
        }
        return export_function_dict[self.data_export_request.data_format.value]

    def _make_zip_file(self) -> io.BytesIO:
        """
        Creates a zip file containing all blobs stored in the specified container in Azure Blob Storage.

        Returns:
            io.BytesIO: A byte stream representing the created zip file.
        """
        # Get the container with all the export files
        container_client = (
            services.file_attachments.create_or_get_existing_container_client(
                self.blob_service_client, self.container_name
            )
        )
        blob_list = container_client.list_blobs()
        # now we need to make a zip file from these files
        zip_stream = io.BytesIO()
        with zipfile.ZipFile(zip_stream, "w", zipfile.ZIP_DEFLATED) as zip_object:
            for blob in blob_list:
                # Adding files that need to be zipped
                fileblob = container_client.download_blob(blob.name)
                data = fileblob.readall()
                # Adding the data to the zip file
                zip_object.writestr(blob.name, data)
        return zip_stream

    def _update_requestor(self):
        """Sends a notification and an email to the requestor telling them their export is ready."""
        models.Notification.create(
            person=self.data_export_request.requestor_person.id,
            type="Data Export Ready",
            info="Your data export is ready.<br>"
            "To respond please follow the link sent to your email account. This request will expire after 7 days.",
            time=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
            status="active",
        )
        services.email_services.send_data_export_ready_message(
            self.data_export_request.requestor_person.user, self.data_export_request
        )

    def _update_db_with_hash(self, checksum: str):
        """
        Adds the checksum for an export file to the database
        Args:
            checksum - the md5 checksum for an exported file.
        """
        self.data_export_request.hash = checksum
        db.session.commit()

    def _upload_zip(self, zip_stream: io.BytesIO) -> BlobClient:
        """
        Uploads a zip file stream to Azure Blob Storage and returns the BlobClient for the uploaded blob.
        Blobs in this container should have a time-based deletion rule applied to them through the Azure storage account

        Args:
            zip_stream (io.BytesIO): The byte stream of the zip file to be uploaded.

        Returns:
            BlobClient: A BlobClient instance representing the uploaded blob.
        """
        # we name the zip file with the uuid, to ensure it has a unique name within the export container.
        file_ext = get_export_file_extension(self.data_export_request.data_format.value)
        zip_blob_name = self.data_export_request.uuid + file_ext
        blob_client = services.file_attachments.get_blob_client(
            "export-outputs", zip_blob_name
        )
        upload_blob = io.BytesIO(zip_stream.getvalue())
        blob_client.upload_blob(upload_blob)
        return blob_client

    @staticmethod
    def _confirm_upload(blob_client: BlobClient, original_hash: str):
        """
        Confirms the blob 1) exists. 2) the checksums match when downloaded

        Args:
            blob_client - the client for the zip blob we are testing has uploaded correctly
            original_hash - the hash generated when we made the zip file

        """
        # confirm upload exists
        if not blob_client.exists():
            return False
        # confirm checksum
        download_stream = blob_client.download_blob()
        checksum = hashlib.md5(download_stream.readall()).hexdigest()
        if not checksum == original_hash:
            return False
        return True

    def _make_eln_export(self):
        eln_export = services.data_export.eln_file.ELNFileExport(
            self.data_export_request,
        )
        eln_export.make_eln_file()
        self._make_zip()

    def _make_rdf_export(self):
        """Creates a zip file of reactions saved as RDFs or JSONs if no SMILES are present"""
        # iterate through the reactions and save each as an azure blob
        for reaction in self.data_export_request.reactions:
            # save all files to the container
            rdf = services.data_export.serialisation.ReactionDataFileExport(
                reaction, reaction.reaction_id, self.data_export_request.uuid
            )
            rdf.save()
        self._make_zip()

    def _make_rxn_export(self):
        """Creates a zip file of reactions saved as RXNs if SMILES are present"""
        # iterate through the reactions and save each as an azure blob
        for reaction in self.data_export_request.reactions:
            # save all files to the container
            rxn = services.data_export.serialisation.RXNFileExport(
                reaction, reaction.reaction_id, self.data_export_request.uuid
            )
            rxn.save()
        self._make_zip()

    def _make_json_export(self):
        """Creates a zip file of reactions saved as JSONs"""
        # iterate through the reactions and save each as an azure blob
        for reaction in self.data_export_request.reactions:
            # save all files to the container
            json_file = services.data_export.serialisation.JsonExport(
                reaction, reaction.reaction_id, self.data_export_request.uuid
            )
            json_file.save()
        self._make_zip()

    def _make_csv_export(self):
        """Create a single CSV to describe all reactions in export"""
        df_rows = []
        csv_export = services.data_export.tabular.CsvExport(self.data_export_request)
        for reaction in self.data_export_request.reactions:
            df_rows.append(csv_export.make_row(reaction))
        # save export df as blob as .surf
        df = pd.concat(df_rows, ignore_index=True)
        csv_contents = services.data_export.tabular.make_csv(df)
        # save the csv file with all the export reactions to Azure.
        file_ext = get_export_file_extension(self.data_export_request.data_format.value)
        save_blob(
            "export-outputs", self.data_export_request.uuid + file_ext, csv_contents
        )

    def _make_surf_export(self):
        """
        Create a single csv to describe all reactions in export with standardised headings
        """
        df_rows = []
        surf_export = services.data_export.tabular.SurfExport(self.data_export_request)
        for reaction in self.data_export_request.reactions:
            df_rows.append(surf_export.make_row(reaction))
        df = pd.concat(df_rows, ignore_index=True)
        csv_contents = services.data_export.tabular.make_csv(df)
        # save the surf file with all the export reactions to Azure.
        file_ext = get_export_file_extension(self.data_export_request.data_format.value)
        save_blob(
            "export-outputs", self.data_export_request.uuid + file_ext, csv_contents
        )

    def _make_pdf_export(self):
        """
        Creates a zip file of reactions saved as PDFs
        """
        # the PDFs are autogenerated when making a reaction summary.
        # we need to iterate through and copy each blob to the export container
        blob_service_client = services.file_attachments.connect_to_azure_blob_service()
        services.file_attachments.create_or_get_existing_container_client(
            blob_service_client, self.data_export_request.uuid
        )
        for reaction in self.data_export_request.reactions:
            self._copy_autogenerated_pdf(reaction, blob_service_client)
        self._make_zip()

    def _copy_autogenerated_pdf(
        self, reaction: models.Reaction, blob_service_client: BlobServiceClient
    ):
        """
        Copies the autogenerated pdf of a reaction, if it exists, to the export container
        Args:
            reaction - the reaction we are trying to copy the autogenerated PDF of
            blob_service_client - the blob service client which is connected to Azure, for all containers.

        """
        # get the blob
        autogenerated_pdf_blob_name = (
            services.file_attachments.get_autogenerated_blob_uuid(reaction.id)
        )
        if autogenerated_pdf_blob_name:
            source_blob_client = blob_service_client.get_blob_client(
                container=current_app.config["STORAGE_CONTAINER"],
                blob=autogenerated_pdf_blob_name,
            )
            # Get a BlobClient for the destination blob
            new_blob = blob_service_client.get_blob_client(
                self.data_export_request.uuid, reaction.reaction_id + ".pdf"
            )
            # Copy the blob from source to destination
            new_blob.start_copy_from_url(source_blob_client.url)


def save_blob(container_name: str, filename: str, file_contents: bytearray):
    """
    Saves the blob to Azure blob service
    Args:
        container_name - the container the blob is being saved in
        filename - the name of the blobfile
        file_contents - the file data

    """
    blob_client = services.file_attachments.get_blob_client(container_name, filename)
    # Upload the blob data
    upload = io.BytesIO(file_contents)
    blob_client.upload_blob(upload, blob_type="BlockBlob")
    # confirm upload
    if not blob_client.exists():
        print(f"blob {filename} upload failed")
        abort(401)
