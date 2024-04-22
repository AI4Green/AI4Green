import hashlib
import io
import zipfile
from datetime import datetime, timedelta
from typing import Callable

import pytz
from azure.storage.blob import BlobClient, BlobSasPermissions, generate_blob_sas
from flask import current_app
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
    blob_client = blob_service_client.get_blob_client(
        container="export-outputs", blob=data_export_request.uuid + ".zip"
    )
    # Set the start time and expiry time for the SAS token
    start_time = datetime.utcnow()
    expiry_time = datetime.utcnow() + timedelta(minutes=2)

    # Define permissions for the SAS token
    permissions = BlobSasPermissions(read=True)

    # Generate the SAS token
    sas_token = generate_blob_sas(
        blob_client.account_name,
        blob_client.container_name,
        blob_client.blob_name,
        account_key=blob_service_client.credential.account_key,
        permission=permissions,
        expiry=expiry_time,
        start=start_time,
    )

    # Construct the SAS URL
    return f"{blob_client.url}?{sas_token}"


def initiate(current_app_context: AppContext, data_export_request_id: int):
    """
    Calls the DataExport class to start making the export file.
    Called from here to instantiate the class object during the threaded process.
    """
    with current_app_context:
        export = DataExport(data_export_request_id)
        export.create()


class DataExport:
    """Class to control the creation the data export, once approved."""

    # def __init__(self, data_export_request: models.DataExportRequest):
    def __init__(self, data_export_request_id: int):
        # db.session.refresh(data_export_request)
        self.data_export_request = models.DataExportRequest.query.get(
            data_export_request_id
        )
        self.container_name = self.data_export_request.uuid
        self.blob_service_client = (
            services.file_attachments.connect_to_azure_blob_service()
        )

    def create(self):
        """
        Main function where we decide the appropriate export type to make and then create a zip file
        and notify the requestor once this is complete.
        """
        export_function = self._get_export_function()
        export_function()
        # Make the zip file and generate an md5 checksum to check file integrity
        zip_stream = self._make_zip_file()
        checksum = hashlib.md5(zip_stream.getvalue()).hexdigest()
        # upload file and confirm the upload was successful
        blob_client = self._upload_zip(zip_stream)
        self._confirm_upload(blob_client, checksum)
        # update database entry with hash
        self._update_db_with_hash(checksum)
        self._update_requestor()
        # delete container
        self.blob_service_client.delete_container(self.container_name)
        # the zip file will be deleted after 7 days. This is controlled by Azure Blob Storage Lifecycle policies.

    def _get_export_function(self) -> Callable:
        """Returns the function which matches the requested data export format"""
        export_function_dict = {
            "RDF": self._make_rdf_export,
            "CSV": self._make_csv_export,
            "JSON": self._make_json_export,
            "SURF": self._make_surf_export,
            "PDF": self._make_pdf_export,
            "ELN": self._make_eln_export,
        }
        return export_function_dict[self.data_export_request.data_format.value]

    def _make_zip_file(self):
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
        services.email.send_data_export_ready_message(
            self.data_export_request.requestor_person.user, self.data_export_request
        )

    def _update_db_with_hash(self, checksum: str):
        """Adds the checksum for an export file to the database"""
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
        zip_blob_name = self.data_export_request.uuid + ".zip"
        blob_client = services.file_attachments.get_blob_client(
            "export-outputs", zip_blob_name
        )
        upload_blob = io.BytesIO(zip_stream.getvalue())
        blob_client.upload_blob(upload_blob, overwrite=True)
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

    def _make_rdf_export(self):
        """
        The main function called to create a zip file of reactions saved as RDFs or JSONs if no SMILES are present
        """
        # iterate through the reactions and save each as an azure blob
        for reaction in self.data_export_request.reactions:
            # save all files to the container
            rdf = services.data_export.reaction_data_file.ReactionDataFile(
                reaction, reaction.reaction_id, self.data_export_request.uuid
            )
            rdf.save()

    def _make_csv_export(self):
        pass

    def _make_json_export(self):
        pass

    def _make_surf_export(self):
        pass

    def _make_pdf_export(self):
        """
        The main function called to create a zip file of reactions saved as PDFs
        """
        # the PDFs are autogenerated when making a reaction summary. - may change in future for most up-to date results.
        # we need to iterate through and copy each blob to the container for this export
        blob_service_client = services.file_attachments.connect_to_azure_blob_service()
        services.file_attachments.create_or_get_existing_container_client(
            blob_service_client, self.data_export_request.uuid
        )
        for reaction in self.data_export_request.reactions:
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
                    self.data_export_request.uuid, reaction.reaction_id
                )
                # Copy the blob from source to destination
                new_blob.start_copy_from_url(source_blob_client.url)

    def _make_eln_export(self):
        pass
