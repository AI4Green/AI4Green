from dataclasses import asdict, dataclass
from datetime import datetime

from flask import current_app, json
from sources.services.message_queue import MessageSerialiserMixin
from sources import services


@dataclass
class ReactionEditMessage(MessageSerialiserMixin):
    """Class for creating a kafka message for the reaction_editing_history topic"""

    person: int
    workgroup: int
    workbook: int
    reaction: int
    field_name: str
    change_details: dict
    date: str

    def serialise(self):
        """Convert a message into a JSON object with a schema and payload.

        The schema is an object that lists the fields and their types.
        The payload is an object representing the message class in JSON format.

        Returns:
            str: The message class formated with shcema and payload.
        """
        schema = {
            "type": "struct",
            "optional": False,
            "fields": [
                {"field": "person", "type": "int32"},
                {"field": "workgroup", "type": "int32"},
                {"field": "workbook", "type": "int32"},
                {"field": "reaction", "type": "int32"},
                {"field": "field_name", "type": "string"},
                {
                    "field": "change_details",
                    "type": "map",
                    "keys": "string",
                    "values": "string",
                },
                {"field": "date", "type": "string"},
            ],
        }
        payload = asdict(self)
        serialised = {"schema": schema, "payload": payload}
        return serialised


def send_message(message: ReactionEditMessage):
    """Send a message to the kafka producer in the reaction_editing_history topic.
    Args:
        message (ReactionEditMessage): The message to send to the queue in the ReactionEditMessage format
    """
    producer = current_app.config["MESSAGE_QUEUE_PRODUCER"]
    producer.send("reaction_editing_history", message.serialise())


def add_new_reaction(person, workbook, reaction_id, reaction_name):
    """Record that a new reaction was created"""
    # get reaction.id from reaction.reaction_id
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook.id
    )
    change_details = {"reaction_name": reaction_name}

    message = ReactionEditMessage(
        person=person.id,
        workgroup=workbook.group,
        workbook=workbook.id,
        reaction=reaction.id,
        field_name="New Reaction",
        change_details=change_details,
        date=datetime.now().strftime("%Y-%m-%d"),
    )
    send_message(message)


def edit_reaction(person, reaction, old_reaction_details, new_reaction_details):
    """Record that a reaction was edited through the sketcher autosave or other means"""
    diff = get_differences(old_reaction_details, new_reaction_details)
    if diff:
        message = ReactionEditMessage(
            person=person.id,
            workgroup=reaction.workbook.group,
            workbook=reaction.workbook.id,
            reaction=reaction.id,
            field_name="Edited Reaction",
            change_details=diff,
            date=datetime.now().strftime("%Y-%m-%d"),
        )
        send_message(message)


def autosave_reaction(person, reaction, old_reaction_details):
    """Record that a reaction was edited through autosave"""
    new_reaction_details = services.reaction.get_reaction_details(reaction)

    # convert json fields to dicts
    old_data_normalised = normalise_json_fields(old_reaction_details)
    new_data_normalised = normalise_json_fields(new_reaction_details)

    # find difference
    diff = get_differences(old_data_normalised, new_data_normalised, exclude_paths=True)
    if diff:
        message = ReactionEditMessage(
            person=person.id,
            workgroup=reaction.workbook.group,
            workbook=reaction.workbook.id,
            reaction=reaction.id,
            field_name="Edited Reaction",
            change_details=diff,
            date=datetime.now().strftime("%Y-%m-%d"),
        )
        send_message(message)


def clone_reaction(person, workbook, new_reaction, old_reaction):
    """Record that a reaction was cloned"""
    change_details = {
        "reaction_id": {
            "old_value": old_reaction.reaction_id,
            "new_value": new_reaction.reaction_id,
        },
        "reaction_name": {"old_value": None, "new_value": new_reaction.name},
    }

    message = ReactionEditMessage(
        person=person.id,
        workgroup=workbook.group,
        workbook=workbook.id,
        reaction=old_reaction.id,
        field_name="Cloned Reaction",
        change_details=change_details,
        date=datetime.now().strftime("%Y-%m-%d"),
    )
    send_message(message)


def delete_reaction(reaction):
    """Record that a reaction was deleted"""
    message = ReactionEditMessage(
        person=reaction.creator_person.id,
        workgroup=reaction.workbook.group,
        workbook=reaction.workbook.id,
        reaction=reaction.id,
        field_name="Deleted Reaction",
        change_details={},
        date=datetime.now().strftime("%Y-%m-%d"),
    )
    send_message(message)


def upload_file(person, reaction, file_names):
    """Record that one or more experimental files were added to a reaction"""
    change_details = {"file_names": file_names}
    message = ReactionEditMessage(
        person=person.id,
        workgroup=reaction.workbook.group,
        workbook=reaction.workbook.id,
        reaction=reaction.id,
        field_name="Uploaded Files",
        change_details=change_details,
        date=datetime.now().strftime("%Y-%m-%d"),
    )
    send_message(message)


def delete_file(person, reaction, file_name):
    """Record that an experimental file was removed from a reaction"""
    change_details = {"file_names": file_name}
    message = ReactionEditMessage(
        person=person.id,
        workgroup=reaction.workbook.group,
        workbook=reaction.workbook.id,
        reaction=reaction.id,
        field_name="Deleted Files",
        change_details=change_details,
        date=datetime.now().strftime("%Y-%m-%d"),
    )
    send_message(message)


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
    """Identify the differences between two dicts, and return with nested keys 'old_value' and 'new_value'.
    Args:
        old_dict (dict): original dict
        new_dict (dict): updated dict
        exclude_paths (bool, optional): whether to exclude all paths not listed in included_paths()
        parent_key (str, optional): only used internally, for recursive calls to nested keys
    """
    diffs = {}

    all_keys = set(old_dict.keys()) | set(new_dict.keys())

    for key in all_keys:
        if exclude_paths and key not in included_paths():
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
