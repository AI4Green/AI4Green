from flask import json
from sources import models, services
from sources.extensions import db


def add(person, workbook, reaction_id, field_name, change_details):
    """
    Add an entry to the ReactionEditingHistory table to record a change in a reaction.
    """

    models.ReactionEditingHistory.create(
        person_id=person.id,
        workbook_id=workbook.id,
        reaction_id=reaction_id,
        field_name=field_name,
        change_details=change_details,
    )


def get_history_from_reaction(reaction_id: int):
    """
    Returns the reaction editing history for reaction
    """
    return (
        db.session.query(models.ReactionEditingHistory)
        .filter(models.ReactionEditingHistory.reaction_id == reaction_id)
        .all()
    )


def add_new_reaction(person, workbook, reaction_id, reaction_name):
    """Record that a new reaction was created"""
    # get reaction.id from reaction.reaction_id
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook.id
    )

    field_name = "New Reaction"
    change_details = {"reaction_name": reaction_name}  # TODO: come back to this
    add(person, workbook, reaction.id, field_name, json.dumps(change_details))


def edit_reaction(person, reaction, old_reaction_details, new_reaction_details):
    """Record that a reaction was edited through the sketcher autosave or other means"""
    field_name = "Edited Reaction"
    diff = get_differences(old_reaction_details, new_reaction_details)
    if diff:
        change_details = json.dumps(diff)
        add(person, reaction.workbook, reaction.id, field_name, change_details)


def autosave_reaction(person, reaction, old_reaction_details):
    """Record that a reaction was edited through autosave"""
    new_reaction_details = services.reaction.get_reaction_details(reaction)

    # convert json fields to dicts
    old_data_normalised = normalise_json_fields(old_reaction_details)
    new_data_normalised = normalise_json_fields(new_reaction_details)

    # find difference
    diff = get_differences(old_data_normalised, new_data_normalised, exclude_paths=True)
    if diff:
        field_name = "Edited Reaction"
        change_details = json.dumps(diff)

        add(person, reaction.workbook, reaction.id, field_name, change_details)


def clone_reaction(person, workbook, new_reaction, old_reaction):
    """Record that a reaction was cloned"""
    field_name = "Cloned Reaction"
    change_details = {
        "reaction_id": {
            "old_value": old_reaction.reaction_id,
            "new_value": new_reaction.reaction_id,
        },
        "reaction_name": {"old_value": None, "new_value": new_reaction.name},
    }
    add(person, workbook, old_reaction.id, field_name, json.dumps(change_details))


def delete_reaction(reaction):
    """Record that a reaction was deleted"""
    add(reaction.creator_person, reaction.workbook, reaction.id, "Deleted Reaction", "")


def upload_file(person, reaction, file_names):
    """Record that one or more experimental files were added to a reaction"""
    field_name = "Uploaded Files"
    change_details = {"file_names": file_names}
    add(person, reaction.workbook, reaction.id, field_name, json.dumps(change_details))


def delete_file(person, reaction, file_name):
    """Record that an experimental file was removed from a reaction"""
    field_name = "Deleted File"
    change_details = {"file_name": file_name}
    add(person, reaction.workbook, reaction.id, field_name, json.dumps(change_details))


def normalise_json_fields(data):
    """converts dict with json strings to nested dict"""
    for key, value in data.items():
        if isinstance(value, str):
            try:
                json_value = json.loads(value)
                data[key] = json_value
            except (ValueError, TypeError):
                pass
    return data


def get_differences(old_dict, new_dict, exclude_paths=False, parent_key=""):
    diffs = {}

    all_keys = set(old_dict.keys()) | set(new_dict.keys())

    for key in all_keys:
        if exclude_paths and key not in included_paths():
            print(key)
            continue
        full_key = f"{parent_key}.{key}" if parent_key else key
        old_value = old_dict.get(key, None)
        new_value = new_dict.get(key, None)

        if isinstance(old_value, dict) and isinstance(new_value, dict):
            nested_diffs = get_differences(old_value, new_value, full_key)
            diffs.update(nested_diffs)
        elif old_value != new_value:
            diffs[full_key] = {"old_value": old_value, "new_value": new_value}

    return diffs


def included_paths():
    # exclude auto-generated data
    paths = [
        "complete",
        "reaction_smiles",
        "description",
        "reactants",
        "products",
        "reagents",
        "solvent",
        "polymerisation_type",
        "reaction_table_data",
        "summary_table_data",
        # reaction_table_data:
        "reactant_names",
        "reactant_masses",
        "reactant_equivalents",
        "reactant_physical_forms_text",
        "limiting_reactant_table_number",
        "reagent_names",
        "reagent_equivalents",
        "reagent_physical_forms_text",
        "reagent_masses",
        "solvent_names",
        "solvent_volumes",
        "solvent_physical_forms_text",
        "product_names",
        "product_equivalents",
        "product_masses",
        "product_physical_forms_text",
        "main_product",
        "amount_units",
        "mass_units",
        "volume_units",
        "solvent_volume_units",
        "product_amount_units",
        "product_mass_units",
        "reactant_mns",
        "product_mns",
        # summary_table_data:
        "real_product_mass",
        "unreacted_reactant_mass",
        "polymer_mn",
        "polymer_mw",
        "polymer_dispersity",
        "polymer_mass_method",
        "polymer_mass_calibration",
        "polymer_tg",
        "polymer_tm",
        "polymer_tc",
        "polymer_thermal_method",
        "polymer_thermal_calibration",
        "reaction_temperature",
        "batch_flow",
        "element_sustainability",
        "isolation_method",
        "catalyst_used",
        "catalyst_recovered",
        "custom_protocol1",
        "custom_protocol2",
        "other_hazards_text",
        "researcher",
        "supervisor",
        "radio_buttons",
    ]
    return paths
