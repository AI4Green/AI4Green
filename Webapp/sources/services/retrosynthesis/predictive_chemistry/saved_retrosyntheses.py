import json
from typing import Dict, List, Tuple

import dash_bootstrap_components as dbc
from dash import callback_context, ctx, html
from flask_login import current_user
from sources import models, services
from sources.extensions import db
from sqlalchemy import func

from .utils import rdkit_smiles_to_image


class SaveRetrosynthesis:
    """Validates and saves a retrosynthesis to the database. Updates the tracker number when finished"""

    def __init__(
        self,
        name: str,
        solved_routes: dict,
        conditions: dict,
        sustainability: dict,
        workbook_id: int,
        new_retrosynthesis_saved_tracker: int,
        retrosynthesis_uuid: str,
    ):
        """
        Creates a new instance of the SaveRetrosynthesis class
        Args:
             name - the name of the retrosynthesis the user sees and in the database
             solved_routes - the dictionary of retrosynthetic routes
             conditions - the dictionary of conditions data
             sustainability - the dictionary of sustainability data
             workbook_id - the primary key ID of the workbook the retrosynthetisis is being saved in
             new_retrosynthesis_saved_tracker - the number which tracks the number of saves made in this current session
             retrosynthesis_uuid - a unique identifier for the retrosynthesis
        """
        self.name = name
        self.solved_routes = solved_routes
        self.conditions = conditions
        self.sustainability = sustainability
        self.workbook_id = workbook_id
        self.new_retrosynthesis_saved_tracker = new_retrosynthesis_saved_tracker
        self.retrosynthesis_uuid = retrosynthesis_uuid
        self.user_message = ""
        self.validation = ""

    def save_process(self) -> Tuple[str, int]:
        """
        Validates the save and if successful calls functions to save to database and returns message to user
        Returns:
            user message as a string
            updated save tracker
        """
        user_message, validation = self.validate_save()
        if validation == "success":
            self.save_retrosynthesis_to_db()
            self.new_retrosynthesis_saved_tracker += 1
        return user_message, self.new_retrosynthesis_saved_tracker

    def validate_save(self) -> Tuple[str, str]:
        """
        Validates the name and checks for it being a non-duplicate
        Returns:
            string either 'success' or 'failed'
            a message for the user with more information
        """
        name_validation = self.validate_name()

        uuid_validation, duplicate_retrosynthesis_name = self.validate_non_duplicate()
        if name_validation == "failed" or uuid_validation == "failed":
            validation = "failed"
            if name_validation == "failed":
                user_message = (
                    "A retrosynthesis with this name already exists in this workbook. "
                    "Please use a different name."
                )
            if uuid_validation == "failed":
                user_message = f"This route has already been saved with the name {duplicate_retrosynthesis_name}"
        else:
            user_message = (f"Retrosynthesis: {self.name} saved successfully",)
            validation = "success"

        return user_message, validation

    def validate_name(self) -> str:
        """
        Checks the save name is unique within the workbook

        Returns:
            'success' or 'failed'
        """
        unique_name_check = (
            db.session.query(models.Retrosynthesis)
            .filter(func.lower(models.Retrosynthesis.name) == self.name.lower())
            .join(models.WorkBook)
            .filter(models.WorkBook.id == self.workbook_id)
            .first()
        )
        if unique_name_check:
            return "failed"
        return "success"

    def validate_non_duplicate(self) -> Tuple[str, str]:
        """
        Checks the uuid is unique. Prevents user saving same retrosynthesis twice in the same workbook.
        Returns:
            'success' or 'failed'
            if failed, also returns the name

        """
        unique_check = (
            db.session.query(models.Retrosynthesis)
            .filter(models.Retrosynthesis.uuid == self.retrosynthesis_uuid)
            .join(models.WorkBook)
            .filter(models.WorkBook.id == self.workbook_id)
            .first()
        )
        if unique_check:
            return "failed", unique_check.name
        return "success", ""

    def save_retrosynthesis_to_db(self):
        """Converts the JSONS to strings and saves the retrosynthesis to the database"""
        target_smiles = self.solved_routes["routes"]["Route 1"]["steps"][0]["smiles"]
        solved_routes_json = json.dumps({"routes": self.solved_routes})
        conditions_json = json.dumps({"routes": self.conditions})
        sustainability_json = json.dumps({"routes": self.sustainability})
        add(
            self.name,
            target_smiles,
            self.retrosynthesis_uuid,
            self.workbook_id,
            conditions_json,
            sustainability_json,
            solved_routes_json,
        )


def make_retrosynthesis_card_list(selected_workbook_id: int) -> html.Div:
    """
    Makes a HTML card for each saved retrosynthesis records within the chosen workbook
    Args:
        selected_workbook_id - the workbook we are checking the user has access to
    Returns:
        a html div containing all of the cards
    """
    retrosynthesis_list = list_from_workbook(selected_workbook_id)
    card_list = []
    for idx, retrosynthesis in enumerate(retrosynthesis_list):
        img_data = rdkit_smiles_to_image(retrosynthesis.target_smiles)
        card_list.append(
            dbc.Card(
                className="mb-4 card-body",
                children=[
                    html.Div(
                        className="pl-3 pt-1 pb-1",
                        style={"margin-bottom": "-1rem"},
                        children=[
                            html.H4(retrosynthesis.name),
                            html.Div(
                                children=[
                                    html.P(retrosynthesis.creator_person.user.fullname),
                                    html.P(
                                        str(retrosynthesis.time_of_creation)[:-7],
                                        className="small text-muted",
                                    ),
                                ]
                            ),
                        ],
                    ),
                    html.Img(
                        src=img_data,
                        style={"background-color": "transparent", "opacity": "100"},
                    ),
                    html.Button(
                        "Reload",
                        className="btn-primary",
                        value=retrosynthesis.id,
                        n_clicks=0,
                        id={"type": "retrosynthesis-reload", "index": idx},
                    ),
                ],
            )
        )
    card_group = html.Div(
        children=card_list,
        className="card shadow-0 border",
        style={"margin-left": "1rem"},
    )
    return card_group


def get_retrosynthesis_to_reload_id(reload_id_values: List[int]) -> int:
    """
    Takes the list of saved retrosynthesis ID the user has access to and uses the context to get the ID of the one
    they clicked using the index

    Args:
        reload_id_values - list of all saved retrosynthesis IDs in the selected workbook in the dropdown
    Returns:
        retrosynthesis_to_reload_id - the database ID of the retrosynthesis to reload.

    """
    retrosynthesis_to_reload = ctx.triggered_id
    retrosynthesis_to_reload_id = None
    if retrosynthesis_to_reload:
        index = retrosynthesis_to_reload["index"]
        retrosynthesis_to_reload_id = reload_id_values[index]
    return retrosynthesis_to_reload_id


def assert_button_clicked(reload_button_clicks: List) -> bool:
    """
    Checks whether the save button has been clicked.
    Args:
        reload_button_clicks - a list with values for the different button clicks if all 0, means unclicked.
    Returns:
        True if save button was clicked.
    """
    # check button has been clicked to prevent firing on initial load
    zero_clicks = all(v == 0 for v in reload_button_clicks)
    changed_ids = [p["prop_id"] for p in callback_context.triggered][0]
    if "n_clicks" in changed_ids and not zero_clicks:
        return True
    return False


def get_reloaded_retrosynthesis(
    retrosynthesis_to_reload_id: int,
) -> Tuple[Dict, Dict, Dict, str]:
    """
    When user clicks reload on a saved retrosynthesis, this function returns the relevant data as JSON
    Args:
        retrosynthesis_to_reload_id - the primary key integer id of the retrosynthesis being relaode.d
    Returns:
        The retrosynthetic routes
        The conditions data
        The sustainability data
        The unique identifier for this retrosynthesis
    """
    retrosynthesis_to_reload = get(retrosynthesis_to_reload_id)
    routes = json.loads(retrosynthesis_to_reload.routes)["routes"]
    # route_dict_array = str_to_dict_array(routes)
    conditions = json.loads(retrosynthesis_to_reload.conditions)["routes"]
    # condition_dict_array = str_to_dict_array(conditions)
    sustainability = json.loads(retrosynthesis_to_reload.sustainability)["routes"]
    retrosynthesis_uuid = retrosynthesis_to_reload.uuid
    return routes, conditions, sustainability, retrosynthesis_uuid


def save_new_reaction_from_retrosynthesis(
    workbook: models.WorkBook,
    reaction_name: str,
    reaction_id: str,
    reaction_smiles: str,
) -> str:
    """
    Makes a new reaction from a retrosynthetic reaction.
    Args:
        workbook - the workbook the reaction is being made in
        reaction_name - the user-provided name of the reaction being made
        reaction_id - the autogenerated id such as DW1-001
        reaction_smiles - the SMILES for the reaction
    """
    # finds workgroup object (needs institution later)
    name_check = check_reaction(workbook, reaction_name)

    creator = current_user.Person

    # check for reaction id - catches errors caused if user has 2 tabs open
    reaction_id_check = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook.id
    )
    if reaction_id_check:
        return "A reaction with this ID already exists. Please refresh the page and try again."

    # if the name check is passed then proceed with making the new reaction
    if "This reaction name is unique" in name_check:
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

        services.reaction.add(
            reaction_name,
            reaction_id,
            creator,
            workbook.id,
            reaction_smiles,
            reaction_table,
            summary_table,
        )
        # load sketcher
        return "New reaction made"
    else:
        return name_check


def check_reaction(workbook: models.WorkBook, reaction_name: str) -> str:
    """
    Checks the reaction name is valid and unique.
    Args:
        workbook - the workbook the reaction is being created in
        reaction_name - the user-provided name of the reaction being made

    """
    # Tells the user they must give the reaction a name to save it
    if not reaction_name:
        return "The reaction must be given a name"
    if not reaction_name.replace(" ", "").replace("-", "").isalnum():
        return "Reaction names cannot contain special characters, only letters and numbers!"
    # Tells the user the reaction name must be unique
    reaction_name_check = services.reaction.get_from_name_and_workbook_id(
        reaction_name, workbook.id
    )
    if reaction_name_check is None:
        feedback = "This reaction name is unique"  # added to the reaction database'
    else:
        feedback = "This reaction name is already used. Please choose another name."
    return feedback


def get(retrosynthesis_id: int) -> models.Retrosynthesis:
    """
    Returns a retrosynthesis object using the retrosynthesis id
    Args:
        retrosynthesis_id - the primary key integer id of the retrosynthesis
    Returns:
        the matching retrosynthesis record.
    """
    return (
        db.session.query(models.Retrosynthesis)
        .filter(models.Retrosynthesis.id == retrosynthesis_id)
        .first()
    )


def add(
    name: str,
    target_smiles: str,
    uuid: str,
    workbook_id: int,
    conditions: Dict[str, any],
    sustainability: Dict[str, any],
    routes: Dict[str, any],
) -> models.Retrosynthesis:
    """Saves the new retrosynthesis to the database"""
    creator = current_user.person
    retrosynthesis = models.Retrosynthesis(
        name=name,
        creator=creator,
        workbook=workbook_id,
        conditions=conditions,
        sustainability=sustainability,
        routes=routes,
        target_smiles=target_smiles,
        uuid=uuid,
    )
    db.session.add(retrosynthesis)
    db.session.commit()
    return retrosynthesis


def list_from_workbook(selected_workbook_id: int):
    """
    Lists all the retrosynthesis records for a workbook
    Args:
         selected_workbook_id - the primary key integer id of the workbook we are getting retrosynthetic records for.
    """
    return (
        db.session.query(models.Retrosynthesis)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == selected_workbook_id)
        .all()
    )
