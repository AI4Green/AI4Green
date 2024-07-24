"""
This module receives a reaction name from the reaction
table and saves it in the reaction database
"""
import base64
from datetime import datetime

import pytz
from flask import Response, json, jsonify, request
from flask_api import status
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

from . import save_reaction_bp


@save_reaction_bp.route("/new_reaction", methods=["POST", "GET"])
@login_required
def new_reaction() -> Response:
    """Makes a new reaction after user submits modal window"""
    workbook_name = request.form["workbook"]
    workgroup_name = request.form["workgroup"]

    # finds workbook object (needs institution later)
    workbook_object = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )

    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook_object)
    reaction_name = sanitise_user_input(request.form["reactionName"])
    reaction_id = request.form["reactionID"]

    creator = services.person.from_current_user_email()

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
        services.reaction.add(
            reaction_name,
            reaction_id,
            creator,
            workbook_object.id,
            reaction_table,
            summary_table,
        )
        # load sketcher
        feedback = "New reaction made"
        return jsonify({"feedback": feedback})
    else:
        return name_check


@save_reaction_bp.route("/_autosave", methods=["POST"])
@login_required
def autosave() -> Response:
    """autosave when a field changes in the reaction page"""
    reaction_description = str(request.form["reactionDescription"])
    reaction = services.reaction.get_current_from_request()
    reaction_name = reaction.name

    services.auth.edit_reaction(reaction)

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
    selectivity = request.form["selectivity"]
    conversion = request.form["conversion"]
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
            "selectivity": selectivity,
            "conversion": conversion,
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


@save_reaction_bp.route("/clone_reaction", methods=["POST", "GET"])
@login_required
def clone_reaction() -> Response:
    """
    Takes reactions data from previously saved reaction and copies it into a new reaction which is saved to the database
    real_product_mass and unreacted_reactant_mass are removed from the reaction data for users to input upon reaction completion

    Returns:
        Response as JSON, either New Reaction Made or Reaction Name is not Unique

    """

    old_reaction = services.reaction.get_current_from_request_form()

    new_reaction_name = request.form.get("reactionName")
    workbook_name = request.form.get("workbook")
    workgroup_name = request.form.get("workgroup")
    new_reaction_id = request.form.get("newReactionID")

    if None in [new_reaction_name, workbook_name, workgroup_name, new_reaction_id]:
        return status.HTTP_400_BAD_REQUEST

    workbook_object = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook_object)

    creator = services.person.from_current_user_email()

    remove_yield_dict = json.loads(old_reaction.summary_table_data)
    remove_yield_dict.update({"real_product_mass": "", "unreacted_reactant_mass": ""})

    # check for reaction id - catches errors caused if user has 2 tabs open
    reaction_id_check = services.reaction.get_from_reaction_id_and_workbook_id(
        new_reaction_id, workbook_object.id
    )
    if reaction_id_check:
        feedback = "A reaction with this ID already exists. Please refresh the page and try again."
        return jsonify({"feedback": feedback})

    name_check = check_reaction_name()
    if b"This reaction name is unique" in name_check.data:
        services.reaction.add(
            new_reaction_name,
            new_reaction_id,
            creator,
            workbook_object.id,
            old_reaction.reaction_table_data,
            json.dumps(remove_yield_dict),
            old_reaction.reaction_smiles,
        )
        feedback = "New reaction made"
        return jsonify({"feedback": feedback})
    else:
        return name_check


@save_reaction_bp.route("/_autosave_sketcher", methods=["POST"])
@login_required
def autosave_sketcher() -> Response:
    """Autosave function for saving changes to the sketcher only. Only used before reaction table is generated."""

    reaction = services.reaction.get_current_from_request()
    services.auth.edit_reaction(reaction)

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
    services.auth.reaction_files(permission_level="edit")
    new_upload = services.file_attachments.UploadExperimentDataFiles(request)
    new_upload.validate_files()
    new_upload.save_validated_files()
    return jsonify({"uploaded_files": new_upload.uploaded_files})


@save_reaction_bp.route("/view_reaction_attachment", methods=["GET"])
@login_required
def view_reaction_attachment() -> Response:
    """
    Authenticate user has permission to view, then use the uuid to find the file on azure and then return the file.
    """
    services.auth.reaction_files(permission_level="view_only")
    file_uuid = request.args.get("uuid")
    # get the blob file and send the filestream, display name and mimetype to the frontend to be displayed.
    blob_client = services.file_attachments.get_blob(file_uuid)
    file_attachment = services.file_attachments.download_from_blob_client(
        blob_client, file_uuid
    )
    file_object = services.file_attachments.database_object_from_uuid(file_uuid)

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
    services.auth.reaction_files(permission_level="view_only")
    blob_client = services.file_attachments.get_blob()
    file_uuid = request.form["uuid"]
    file_attachment = services.file_attachments.download_from_blob_client(
        blob_client, file_uuid
    )
    file_object = services.file_attachments.database_object_from_uuid(file_uuid)
    mimetype = file_object.file_details["mimetype"]
    display_name = file_object.display_name
    return jsonify(
        {"stream": file_attachment, "mimetype": mimetype, "name": display_name}
    )


@save_reaction_bp.route("/_delete_reaction_attachment", methods=["DELETE"])
@login_required
def delete_reaction_attachment():
    """Delete file attached to reaction"""
    services.file_attachments.delete_file_attachment(request_source="user")
    return "success"
