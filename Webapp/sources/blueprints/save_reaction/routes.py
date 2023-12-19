"""
This module receives a reaction name from the reaction
table and saves it in the reaction database
"""
import base64
import io
import os
import uuid
from datetime import datetime

import magic
import pytz
from azure.core.exceptions import ResourceExistsError
from azure.storage.blob import BlobClient, BlobServiceClient
from flask import Response, abort, current_app, json, jsonify, request
from flask_login import current_user, login_required
from sources import models, services
from sources.auxiliary import (
    abort_if_user_not_in_workbook,
    get_data,
    get_smiles,
    sanitise_user_input,
)
from sources.extensions import db
from sqlalchemy import func
from werkzeug.utils import secure_filename

from . import save_reaction_bp


@save_reaction_bp.route("/new_reaction", methods=["POST", "GET"])
@login_required
def new_reaction() -> Response:
    """Makes a new reaction after user submits modal window"""
    workbook_name = request.form["workbook"]
    workgroup_name = request.form["workgroup"]

    # finds workgroup object (needs institution later)
    workgroup_selected_obj = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )
    workbook_object = (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.name == workbook_name)
        .filter(models.WorkBook.group == workgroup_selected_obj.id)
        .first()
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook_object)
    reaction_name = sanitise_user_input(request.form["reactionName"])
    reaction_id = request.form["reactionID"]

    creator = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .first()
    )
    current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
    # check for reaction id - catches errors caused if user has 2 tabs open
    reaction_id_check = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .filter(models.Reaction.workbooks == workbook_object.id)
        .first()
    )
    if reaction_id_check:
        feedback = "A reaction with this ID already exists. Please refresh the page and try again."
        return jsonify({"feedback": feedback})

    name_check = check_reaction_name()
    # if the name check is passed then proceed with making the new reaction
    if b"This reaction name is unique" in name_check.data:
        # make the reaction table dict with units set to default values
        reaction_table = json.dumps(
            {
                "amount_units": "mmol",
                "mass_units": "mg",
                "volume_units": "mL",
                "solvent_volume_units": "mL",
                "product_amount_units": "mmol",
                "product_mass_units": "mg",
                "reactant_masses": [],
                "reactant_masses_raw": [],
                "reactant_amounts": [],
                "reactant_amounts_raw": [],
                "reactant_volumes": [],
                "reactant_volumes_raw": [],
                "reactant_equivalents": [],
                "reactant_physical_forms": [],
                "reactant_densities": [],
                "reactant_concentrations": [],
                "reagent_names": [],
                "reagent_molecular_weights": [],
                "reagent_densities": [],
                "reagent_concentrations": [],
                "reagent_amounts": [],
                "reagent_amounts_raw": [],
                "reagent_equivalents": [],
                "reagent_physical_forms": [],
                "reagent_hazards": [],
                "reagent_masses": [],
                "reagent_masses_raw": [],
                "reagent_volumes": [],
                "reagent_volumes_raw": [],
                "solvent_volumes": [],
                "solvent_names": [],
                "solvent_concentrations": [],
                "solvent_hazards": [],
                "solvent_physical_forms": [],
                "product_amounts": [],
                "product_amounts_raw": [],
                "product_masses": [],
                "product_masses_raw": [],
                "product_physical_forms": [],
            }
        )

        summary_table = json.dumps(
            {
                "real_product_mass": "",
                "unreacted_reactant_mass": "",
                "reaction_temperature": "",
                "batch_flow": "-select-",
                "element_sustainability": "undefined",
                "isolation_method": "undefined",
                "catalyst_used": "-select-",
                "catalyst_recovered": "-select-",
                "custom_protocol1": "",
                "custom_protocol2": "",
                "other_hazards_text": "",
                "researcher": "",
                "supervisor": "",
                "radio_buttons": [],
            }
        )

        # add reaction to database
        reaction = models.Reaction(
            creator=creator.id,
            reaction_id=reaction_id,
            time_of_creation=current_time,
            name=reaction_name,
            workbooks=workbook_object.id,
            status="active",
            complete="not complete",
            reaction_table_data=reaction_table,
            summary_table_data=summary_table,
        )
        db.session.add(reaction)
        db.session.commit()
        # load sketcher
        feedback = "New reaction made"
        return jsonify({"feedback": feedback})
    else:
        return name_check


def authenticate_user_to_edit_reaction(
    reaction: models.Reaction, file_attachment=False
):
    """
    In addition to frontend validation, backend validation protects against user edited HTML.
    Validates the active user is the creator and validates the reaction is not locked.
    Aborts process with a 400 error if validation is failed
    """
    # validate user is in workbook
    workbook_persons = reaction.workbook.users
    workbook_users = [x.user for x in workbook_persons]
    if current_user not in workbook_users:
        abort(401)
    # validate the user is the creator
    if reaction.creator_person.user.email != current_user.email:
        abort(401)
    if file_attachment:
        return
    # validate the reaction is not locked, unless it is a file attachment being edited.
    if reaction.complete == "complete":
        abort(401)


@save_reaction_bp.route("/_autosave", methods=["POST"])
@login_required
def autosave() -> Response:
    """autosave when a field changes in the reaction page"""
    reaction_description = str(request.form["reactionDescription"])
    reaction = get_current_reaction()
    reaction_name = reaction.name

    authenticate_user_to_edit_reaction(reaction)

    summary_to_print = str(request.form["summary_to_print"])
    current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
    reaction_smiles = str(request.form["reactionSmiles"])

    # reaction table entries
    # find the table number of the limiting reactant e.g js-reactant1
    limiting_reactant = str(request.form["limitingReactantTableNumber"])
    # reactants data from the ajax post
    reactant_primary_keys = get_data("reactantPrimaryKeys")
    reactant_primary_keys_ls = list(filter(None, reactant_primary_keys))
    reactant_smiles_ls = get_smiles(reactant_primary_keys_ls)
    reactant_masses = get_data("reactantMasses")[:-1]
    reactant_masses_raw = get_data("reactantMassesRaw")[:-1]
    reactant_amounts = get_data("reactantAmounts")[:-1]
    reactant_amounts_raw = get_data("reactantAmountsRaw")[:-1]
    reactant_volumes = get_data("reactantVolumes")[:-1]
    reactant_volumes_raw = get_data("reactantVolumesRaw")[:-1]
    reactant_equivalents = get_data("reactantEquivalents")[:-1]
    reactant_physical_forms = get_data("reactantPhysicalForms")[:-1]
    reactant_densities = get_data("reactantDensities")[:-1]
    reactant_concentrations = get_data("reactantConcentrations")[:-1]
    reactant_names = get_data("reactantNames")[:-1]
    reactant_molecular_weights = get_data("reactantMolecularWeights")[:-1]
    reactant_hazards = get_data("reactantHazards")[:-1]
    reactant_physical_forms_text = get_data("reactantPhysicalFormsText")[:-1]
    # reagents data from the ajax post
    reagent_smiles_ls = get_data("reagentSmiles")
    reagent_names = get_data("reagentNames")[:-1]
    reagent_molecular_weights = get_data("reagentMolecularWeights")
    reagent_densities = get_data("reagentDensities")[:-1]
    reagent_concentrations = get_data("reagentConcentrations")[:-1]
    reagent_amounts = get_data("reagentAmounts")[:-1]
    reagent_amounts_raw = get_data("reagentAmountsRaw")[:-1]
    reagent_masses = get_data("reagentMasses")[:-1]
    reagent_masses_raw = get_data("reagentMassesRaw")[:-1]
    reagent_volumes = get_data("reagentVolumes")[:-1]
    reagent_volumes_raw = get_data("reagentVolumesRaw")[:-1]
    reagent_equivalents = get_data("reagentEquivalents")[:-1]
    reagent_physical_forms = get_data("reagentPhysicalForms")[:-1]
    reagent_hazards = get_data("reagentHazards")[:-1]
    reagent_physical_forms_text = get_data("reagentPhysicalFormsText")
    # solvents data from the ajax post
    solvent_primary_keys = get_data("solventPrimaryKeys")
    solvent_primary_keys_ls = list(filter(None, solvent_primary_keys))
    solvent_names = get_data("solventNames")[:-1]
    solvent_concentrations = get_data("solventConcentrations")[:-1]
    solvent_volumes = get_data("solventVolumes")[:-1]
    solvent_physical_forms = get_data("solventPhysicalForms")[:-1]
    solvent_hazards = get_data("solventHazards")[:-1]
    solvent_physical_forms_text = get_data("solventPhysicalFormsText")
    # products data from the ajax post
    main_product = str(request.form["mainProductTableNumber"])
    product_primary_keys = get_data("productPrimaryKeys")
    product_primary_keys_ls = list(filter(None, product_primary_keys))
    product_smiles_ls = get_smiles(product_primary_keys_ls)
    product_physical_form = get_data("productPhysicalForms")[:-1]
    product_amounts = get_data("productAmounts")[:-1]
    product_amounts_raw = get_data("productAmountsRaw")[:-1]
    product_masses = get_data("productMasses")[:-1]
    product_masses_raw = get_data("productMassesRaw")[:-1]
    product_names = get_data("productNames")[:-1]
    product_molecular_weights = get_data("productMolecularWeights")
    product_hazards = get_data("productHazards")
    product_physical_forms_text = get_data("productPhysicalFormsText")
    amount_units = str(request.form["amountUnits"])
    mass_units = str(request.form["massUnits"])
    volume_units = str(request.form["volumeUnits"])
    solvent_volume_units = str(request.form["solventVolumeUnits"])
    product_amount_units = str(request.form["productAmountUnits"])
    product_mass_units = str(request.form["productMassUnits"])
    reaction_table = json.dumps(
        {
            # reaction data
            "reaction_smiles": reaction_smiles,
            "reaction_name": reaction_name,
            "reaction_description": reaction_description,
            # reactant data
            "reactant_smiles": reactant_smiles_ls,
            "reactant_masses": reactant_masses,
            "reactant_amounts": reactant_amounts,
            "reactant_amounts_raw": reactant_amounts_raw,
            "reactant_volumes": reactant_volumes,
            "reactant_volumes_raw": reactant_volumes_raw,
            "reactant_equivalents": reactant_equivalents,
            "reactant_physical_forms": reactant_physical_forms,
            "reactant_masses_raw": reactant_masses_raw,
            "reactant_densities": reactant_densities,
            "reactant_concentrations": reactant_concentrations,
            "limiting_reactant_table_number": limiting_reactant,
            "reagent_smiles": reagent_smiles_ls,
            "reagent_names": reagent_names,
            "reagent_molecular_weights": reagent_molecular_weights,
            "reagent_densities": reagent_densities,
            "reagent_concentrations": reagent_concentrations,
            "reagent_amounts": reagent_amounts,
            "reagent_amounts_raw": reagent_amounts_raw,
            "reagent_equivalents": reagent_equivalents,
            "reagent_physical_forms": reagent_physical_forms,
            "reagent_hazards": reagent_hazards,
            "reagent_masses": reagent_masses,
            "reagent_masses_raw": reagent_masses_raw,
            "reagent_volumes": reagent_volumes,
            "reagent_volumes_raw": reagent_volumes_raw,
            "solvent_ids": solvent_primary_keys_ls,
            "solvent_volumes": solvent_volumes,
            "solvent_names": solvent_names,
            "solvent_concentrations": solvent_concentrations,
            "solvent_hazards": solvent_hazards,
            "solvent_physical_forms": solvent_physical_forms,
            "product_amounts": product_amounts,
            "product_amounts_raw": product_amounts_raw,
            "product_masses": product_masses,
            "product_masses_raw": product_masses_raw,
            "product_smiles": product_smiles_ls,
            "product_physical_forms": product_physical_form,
            "main_product": main_product,
            "amount_units": amount_units,
            "mass_units": mass_units,
            "volume_units": volume_units,
            "solvent_volume_units": solvent_volume_units,
            "product_amount_units": product_amount_units,
            "product_mass_units": product_mass_units,
            "reactant_names": reactant_names,
            "reactant_molecular_weights": reactant_molecular_weights,
            "reactant_hazards": reactant_hazards,
            "reactant_physical_forms_text": reactant_physical_forms_text,
            "reagent_physical_forms_text": reagent_physical_forms_text,
            "solvent_physical_forms_text": solvent_physical_forms_text,
            "product_names": product_names,
            "product_molecular_weights": product_molecular_weights,
            "product_hazards": product_hazards,
            "product_physical_forms_text": product_physical_forms_text,
        }
    )

    # summary table entries
    # product masses and unreacted reactant masses
    real_product_mass = request.form["realProductMass"]
    unreacted_reactant_mass = request.form["unreactedReactantMass"]
    # sustainability data
    reaction_temperature = request.form["reactionTemperature"]
    element_sustainability = request.form["elementSustainability"]
    batch_flow = request.form["batchFlow"]
    isolation_method = request.form["isolationMethod"]
    catalyst_used = request.form["catalystUsed"]
    catalyst_recovered = request.form["catalystRecovered"]
    # radio buttons for standard protocols, disposal, spillage, and hazard categorisation
    radio_buttons = get_data("selectedRadioButtons")[:-1]
    # other hazards textbox, and custom protocol fields
    custom_protocol1 = request.form["customProtocol1"]
    custom_protocol2 = request.form["customProtocol2"]
    other_hazards_text = request.form["otherHazardTextArea"]
    researcher = request.form["researcher"]
    supervisor = request.form["supervisor"]
    mass_efficiency = request.form["massEfficiency"]
    to_export = request.form["toExport"]

    summary_table = json.dumps(
        {
            "real_product_mass": real_product_mass,
            "unreacted_reactant_mass": unreacted_reactant_mass,
            "reaction_temperature": reaction_temperature,
            "element_sustainability": element_sustainability,
            "batch_flow": batch_flow,
            "isolation_method": isolation_method,
            "catalyst_used": catalyst_used,
            "catalyst_recovered": catalyst_recovered,
            "custom_protocol1": custom_protocol1,
            "custom_protocol2": custom_protocol2,
            "other_hazards_text": other_hazards_text,
            "researcher": researcher,
            "supervisor": supervisor,
            "radio_buttons": radio_buttons,
            "summary_to_print": summary_to_print,
            "mass_efficiency": mass_efficiency,
            "to_export": to_export,
        }
    )

    # value is 'complete' if user is trying to lock reaction.
    complete = request.form["complete"]
    feedback = "Reaction Updated!"
    if complete == "complete":
        # check all mandatory fields are complete
        missing_data_fields = [
            reactant_masses,
            reactant_equivalents,
            reagent_equivalents,
            reagent_names,
            solvent_names,
            solvent_volumes,
            reagent_physical_forms,
            solvent_physical_forms,
            product_physical_form,
        ]
        missing_data = False
        for x in missing_data_fields:
            if "" in x or "0" in x:
                missing_data = True
        if unreacted_reactant_mass == "" or real_product_mass == "":
            return jsonify(
                {
                    "feedback": "Please enter unreacted and product mass to mark reaction as complete!"
                }
            )
        elif missing_data:
            return jsonify(
                {
                    "feedback": "Please fill in all required data to mark reaction as complete!"
                }
            )
        else:
            feedback = "Reaction locked"
    # update the database entry for the reaction
    update_dict = {
        "time_of_update": current_time,
        "complete": complete,
        "reaction_smiles": reaction_smiles,
        "description": reaction_description,
        "reactants": reactant_smiles_ls,
        "products": product_smiles_ls,
        "reagents": reagent_smiles_ls,
        "solvent": solvent_primary_keys_ls,
        "reaction_table_data": reaction_table,
        "summary_table_data": summary_table,
    }
    reaction.update(**update_dict)
    return jsonify({"feedback": feedback})


@save_reaction_bp.route("/_autosave_sketcher", methods=["POST"])
@login_required
def autosave_sketcher() -> Response:
    """Autosave function for saving changes to the sketcher only. Only used before reaction table is generated."""

    reaction = get_current_reaction()
    authenticate_user_to_edit_reaction(reaction)

    current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
    reaction_smiles = str(request.form["reactionSmiles"])

    update_dict = {"time_of_update": current_time, "reaction_smiles": reaction_smiles}
    reaction.update(**update_dict)
    feedback = "Reaction Updated!"
    return jsonify({"feedback": feedback})


@save_reaction_bp.route("/_check_reaction", methods=["POST"])
@login_required
def check_reaction_name() -> Response:
    """Checks the reaction name is unique"""
    reaction_name = sanitise_user_input(request.form["reactionName"])
    # Tells the user they must give the reaction a name to save it
    if reaction_name == "":
        feedback = "The reaction must be given a name"
        return jsonify({"feedback": feedback})

    if not reaction_name.replace(" ", "").replace("-", "").isalnum():
        feedback = "Reaction names cannot contain special characters, only letters and numbers!"
        return jsonify({"feedback": feedback})
    workbook_name = str(request.form["workbook"])
    workgroup_name = str(request.form["workgroup"])

    # Tells the user the reaction name must be unique
    reaction_name_check = (
        db.session.query(models.Reaction)
        .filter(func.lower(models.Reaction.name) == reaction_name.lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook_name)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )
    if reaction_name_check is None:
        # populate_database(reaction_name)
        feedback = "This reaction name is unique"  # added to the reaction database'
    else:
        feedback = "This reaction name is already used. Please choose another name."
    return jsonify({"feedback": feedback})


@save_reaction_bp.route("/_upload_experimental_data", methods=["POST"])
@login_required
def upload_experiment_files():
    """Takes a list of files, and saves upon successful validation. Url added to database, file saved to azure blob"""
    authenticate_user(permission_level="edit", request_method="POST")
    new_upload = UploadExperimentDataFiles(request)
    new_upload.validate_files()
    new_upload.save_validated_files()
    return jsonify({"uploaded_files": new_upload.uploaded_files})


@save_reaction_bp.route("/view_reaction_attachment", methods=["GET"])
@login_required
def view_reaction_attachment() -> Response:
    """
    Authenticate user has permission to view, then use the uuid to find the file on azure and then return the file.
    """
    authenticate_user(permission_level="view_only", request_method="GET")
    file_uuid = request.args.get("uuid")
    # get the blob file and send the filestream, display name and mimetype to the frontend to be displayed.
    blob_client = get_blob(file_uuid)
    file_attachment = blob_to_file_attachment(blob_client)
    file_object = file_from_uuid(file_uuid)

    mimetype = file_object.file_details["mimetype"]
    display_name = file_object.display_name
    file_attachment = f"data:{mimetype};base64,{file_attachment}"
    return jsonify(
        {"stream": file_attachment, "mimetype": mimetype, "name": display_name}
    )


@save_reaction_bp.route("/_download_reaction_attachment", methods=["POST"])
@login_required
def download_experiment_files():
    """Take a file and return as attachment to the user"""
    authenticate_user(permission_level="view_only", request_method="POST")
    blob_client = get_blob()
    stream_storage = blob_client.download_blob()
    stream = stream_storage.readall()
    file_attachment = base64.b64encode(stream).decode()
    file_uuid = request.form["uuid"]
    file_object = file_from_uuid(file_uuid)

    mimetype = file_object.file_details["mimetype"]
    display_name = file_object.display_name
    return jsonify(
        {"stream": file_attachment, "mimetype": mimetype, "name": display_name}
    )


@save_reaction_bp.route("/_delete_reaction_attachment", methods=["DELETE"])
@login_required
def delete_reaction_attachment():
    """Delete file attached to reaction"""
    authenticate_user(permission_level="edit", request_method="DELETE")
    # connect to blob
    blob_client = get_blob()
    # delete blob
    blob_client.delete_blob()
    # confirm deletion
    if blob_client.exists():
        abort(401)
    # reflect changes in the database
    file_uuid = request.form["uuid"]
    file_object = file_from_uuid(file_uuid)
    db.session.delete(file_object)
    db.session.commit()
    return "success"


def get_blob(file_uuid: str = None) -> BlobClient:
    # connect to azure
    blob_service_client = connect_to_azure_blob_service()
    # connect to blob
    container_name = current_app.config["STORAGE_CONTAINER"]
    if file_uuid is None:
        file_uuid = request.form["uuid"]
    blob_client = blob_service_client.get_blob_client(
        container=container_name, blob=file_uuid
    )
    return blob_client


def blob_to_file_attachment(blob_client) -> str:
    stream = blob_client.download_blob()
    data = stream.readall()
    file_attachment = base64.b64encode(data).decode()
    return file_attachment


def authenticate_user(permission_level: str, request_method: str):
    """
    Authenticates user to either view or edit the reaction.
    permission_level: Takes value of 'edit' or 'view_only'
    request_method: Value of 'GET' changes behaviour, other strings all have same behaviour (e.g. POST, DELETE)
    """
    if permission_level == "view_only":
        authenticate_user_to_view_files(request_method)
    if permission_level == "edit":
        reaction = get_current_reaction()
        authenticate_user_to_edit_reaction(reaction, file_attachment=True)


def authenticate_user_to_view_files(request_method):
    """Authenticates user as a workbook member or aborts. Gets the workgroup_name, workbook_name, and workbook."""
    if request_method == "GET":
        workgroup_name = request.args.get("workgroup")
        workbook_name = request.args.get("workbook")
    else:
        workgroup_name = request.form["workgroup"]
        workbook_name = request.form["workbook"]
    # validate user belongs to the workbook
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)


def connect_to_azure_blob_service() -> BlobServiceClient:
    """
    Returns the Azure blob service client using the connection string from configs.
    From this object, file attachment contents can be accessed by specifying the container name and file uuid.
    """
    return BlobServiceClient.from_connection_string(
        current_app.config["AZURE_STORAGE_CONNECTION_STRING"]
    )


def file_from_uuid(file_uuid: str) -> models.ReactionDataFile:
    """Return database file object from the uuid. (autogenerated unique identifier)"""
    return (
        db.session.query(models.ReactionDataFile)
        .filter(models.ReactionDataFile.uuid == file_uuid)
        .first()
    )


class UploadExperimentDataFiles:
    def __init__(self, ajax_request):
        self.files_to_upload = ajax_request.files
        self.institution = db.session.query(models.Institution).first()
        self.workgroup = ajax_request.form["workgroup"]
        self.workbook = ajax_request.form["workbook"]
        self.reaction_id = ajax_request.form["reactionID"]
        self.validated_files = []
        self.container_client = None
        self.reaction = get_current_reaction()
        self.uploaded_files = []

    def validate_number_of_attachments(self):
        number_of_attachments = len(self.reaction.file_attachments)
        if number_of_attachments >= 10:
            abort(401)

    def validate_files(self):
        for file in self.files_to_upload.values():
            validation_result = self.file_security_validation(file)
            if validation_result == "success":
                self.validated_files.append(file)
            else:
                abort(401)

    def save_validated_files(self):
        for file in self.validated_files:
            # file name must be unique - needs workgroup/book/reaction_id/filename. Check for uniqueness of filename.
            filename = str(uuid.uuid4())
            self.save_blob(file, filename)
            print("blob uploaded")
            # save file details to database, and shorten filename if too long.
            name, extension = os.path.splitext(file.filename)
            if len(name) > 21:
                file.filename = name[0:10] + "..." + name[-10:-1] + extension
            self.save_file_details_to_database(file, filename)
            self.uploaded_files.append({"uuid": filename, "name": file.filename})

    def save_file_details_to_database(self, file, filename):
        current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        storage_name = self.container_client.account_name
        container_name = self.container_client.container_name
        file_size = os.fstat(file.fileno()).st_size
        file_ext = os.path.splitext(file.filename)[1]
        file_details = {
            "mimetype": file.mimetype,
            "size": file_size,
            "extension": file_ext,
        }
        models.ReactionDataFile.create(
            reaction=self.reaction.id,
            storage_name=storage_name,
            container_name=container_name,
            uuid=filename,
            display_name=file.filename,
            time_of_upload=current_time,
            file_details=file_details,
        )

    def save_blob(self, file, filename):
        # connect to storage account
        blob_service_client = connect_to_azure_blob_service()
        # connect to container
        container_name = current_app.config["STORAGE_CONTAINER"]
        try:
            self.container_client = blob_service_client.create_container(container_name)
        except ResourceExistsError:
            self.container_client = blob_service_client.get_container_client(
                container_name
            )

        blob_client = blob_service_client.get_blob_client(
            container=container_name, blob=filename
        )
        print("uploading blob")
        # Upload the blob data - default blob type is BlockBlob
        file.stream.seek(0)
        upload2 = io.BytesIO(file.stream.read())
        blob_client.upload_blob(upload2, blob_type="BlockBlob")
        # confirm upload
        if not blob_client.exists():
            print(f"blob {filename} upload failed")
            abort(401)

    @staticmethod
    def file_security_validation(file):
        """Validates the filesize, name, extension, and mimetype"""
        filename = file.filename
        secure_filename(filename)
        if filename != "":
            # size under 1 mb
            filesize = os.fstat(file.fileno()).st_size
            if filesize > 1000000:
                print("file too large")
                return "failure"
            # acceptable file extension
            file_ext = os.path.splitext(filename)[1]
            if file_ext not in current_app.config["UPLOAD_EXTENSIONS"]:
                print("incorrect extension")
                return "failure"
            # acceptable mimetype
            mime = magic.Magic(mime=True)
            mime_type = mime.from_buffer(file.read(2048))
            if mime_type not in current_app.config["UPLOAD_MIME_TYPES"]:
                print("incorrect mimetype")
                return "failure"
        else:
            return "failure"
        return "success"


def get_current_reaction():
    reaction_id = str(request.form["reactionID"])
    workgroup_name = str(request.form["workgroup"])
    workbook_name = str(request.form["workbook"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    return (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )
