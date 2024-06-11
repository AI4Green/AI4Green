import json
from typing import List, Optional

from sources import models


def validate_reaction(reaction: models.Reaction) -> bool:
    """Validates a reaction is suitable for export. We use whether the user has entered a yield to decide this.
    Args:
        reaction - the reaction we are validating for export
    Returns:
        True if the reaction is valid (user has entered a yield) and False if it is not.
    """
    yield_mass = json.loads(reaction.summary_table_data).get("real_product_mass")
    if yield_mass and yield_mass.isdigit():
        return True
    return False


class ReactionStringMapper:
    """Collection of static methods to replace strings within the data export methods."""

    @staticmethod
    def css_classes_to_chem21_colour_flag(css_class: str) -> str:
        """Replaces the css classes with their default colour flags - used in elements table"""
        css_chem21_dict = {
            "hazard-reset-hazard": None,
            "hazard-acceptable": "green",
            "hazard-warning": "amber",
            "hazard-hazardous": "red",
        }
        return css_chem21_dict.get(css_class)

    @staticmethod
    def isolation(isolation_method: str) -> Optional[str]:
        """Uses to replace the isolation dropdown index with the corresponding string value"""
        isolation_dict = {
            "0": "",
            "1": "Column",
            "2": "HPLC",
            "3": "Ion exchange",
            "4": "Crystallization",
            "5": "Filtration",
            "6": "Multiple recryst.",
            "7": "Distillation < 140 degC",
            "8": "Distillation > 140 degC",
        }
        return isolation_dict.get(isolation_method)

    @staticmethod
    def element_sustainability(element_value: str) -> str:
        """Replaces the element sustainability dropdown index with the corresponding string value"""
        elements_dict = {
            "3": "+500 years",
            "2": "50-500 years",
            "1": "5-50 years",
            "0": "",
        }
        return elements_dict.get(element_value)

    @staticmethod
    def solvent_sustainability(solvent_flag: int) -> str:
        """Replaces the solvent sustainability numerical flag found in the database with the corresponding string"""
        # get solvent sustainability
        solvent_flag_dict = {
            4: "Recommended",
            3: "Problematic",
            2: "Hazardous",
            1: "Highly Hazardous",
            5: "Unknown",  # non-chem21 elsewhere in code sometimes.
        }
        return solvent_flag_dict.get(solvent_flag)

    @staticmethod
    def physical_forms(physical_form: str):
        """Replaces the physical form dropdown index with the corresponding string value"""
        physical_form_dict = {
            "0": "",
            "1": "Dense Solid",
            "2": "Non-volatile liquid (b.p. > 130 degC)",
            "3": "Unknown",
            "4": "Dusty Solid",
            "5": "Lyophilised solid",
            "6": "Volatile liquid (b.p. 70-130 degC)",
            "7": "Gas",
            "8": "Highly volatile liquid (b.p. < 70 degC)",
            "9": "Aerosol",
            "10": "Solution that promotes skin absorption",
        }
        return physical_form_dict.get(physical_form)

    @staticmethod
    def radio_buttons(radio_button: str) -> str:
        """
        We take the radio button name and change it if it appears in the dict, otherwise we return the original

        Args:
            radio_button - the string from the frontend of the selected radio button

        Returns:
            Either the new string from the dict value or the original if the string was not in the dictionary,
        """
        radio_button_dict = {
            "slight": "hazard: slight",
            "serious": "hazard: serious",
            "major": "hazard: major",
            "low_likelihood": "risk: low likelihood",
            "possible": "risk: possible",
            "frequent_occur": "risk: frequent occur",
            "individual": "consequences: individual",
            "local_labs": "consequences: local labs",
            "building_wide": "consequences: building wide",
        }
        new_string = radio_button_dict.get(radio_button)
        if new_string:
            return new_string
        return radio_button


def get_metadata_keys() -> List[str]:
    """
    Returns the keys used in the metadata
    """
    return [
        "temperature",
        "solvents",
        "reagents",
        "creator_username",
        "creator_workbook",
        "creator_workgroup",
        "time_of_creation",
        "reaction_completed",
        "isolation_method",
        "batch_or_flow",
        "percent_yield",
        "reactants",
        "experimental_writeup",
        "standard_protocols_used",
        "nmr_data",
        "file_attachment_names",
    ]


def remove_default_data(list_: List) -> List:
    """Removes default database values from lists for database reaction components"""
    return [x for x in list_ if x not in ["{", "}", ""]]


def check_compounds_present(db_reaction: models.Reaction, component: str) -> bool:
    """
    Checks if a role of compounds are present within a reaction
    e.g., if component='solvent', a reaction has no solvents returns False
    Args:
        db_reaction - the database entry for the reaction being checked
        component - the role being checked - reactant/reagent/solvent/product
    Returns:
         True if that compound type is present in the reaction.
    """
    return bool(remove_default_data(getattr(db_reaction, component)))
