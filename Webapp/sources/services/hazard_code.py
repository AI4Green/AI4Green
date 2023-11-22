from typing import Dict, List, Tuple

from sources import models
from sources.extensions import db

numerical_rating_from_str = {"": 0, "L": 1, "M": 2, "H": 3, "VH": 4}

str_rating_from_number = {0: "", 1: "L", 2: "M", 3: "H", 4: "VH"}

physical_form_exposure_potential = {  # physical form - risk category dictionary
    "-select-": "M",
    "Dense solid": "L",
    "Non-volatile liquid": "L",
    "Unknown": "M",
    "Dusty Solid": "M",
    "Lyophilised solid": "M",
    "Volatile liquid": "M",
    "Gas": "H",
    "Highly volatile liquid": "H",
    "Aerosol": "H",
    "Solution that promotes skin absorption": "H",
}
hazard_exposure_risk = {  # combination of risks - risk category dictionary
    "VHL": "H",
    "VHM": "VH",
    "VHH": "VH",
    "HL": "M",
    "HM": "H",
    "HH": "H",
    "ML": "L",
    "MM": "M",
    "MH": "H",
    "LL": "L",
    "LM": "L",
    "LH": "M",
}


def get_multiple_compounds_data(
    compounds_hazard_codes_list: List[str], physical_forms_list: List[str]
) -> Tuple[List[int], List[str], List[str], List[str], List[str], List[str], List[str]]:
    """
    For a list of compounds (e.g., products) this returns all their hazard data. Each array has one item per compound

    Args:
        compounds_hazard_codes_list - example: ['H201-H302', 'H401-H310']
        physical_forms_list - example: ['Volatile liquid', 'Dense solid']

    Returns:
        compounds_most_severe_hazard_numerical_ratings - how hazardous each compound is. 4 most hazardous
        compounds_hazard_sentences - each string item is all the h codes and phrases for a compound in a sentence
        compound_hazard_ratings - how hazardous each compound is. VH most hazardous
        compound_hazard_colours - should hazard table cell be white background (hazard-reset) or red (hazard-hazardous)
        compounds_exposure_potentials - the exposure potential of each compound. VH is the highest exposure potential
        compounds_risk_ratings - the risk from a compound. VH is the highest risk
        compounds_risk_colours - should risk table cell be white background(hazard-reset) or red (hazard-hazardous)
    """
    # check data

    # initiate the list variables to be returned
    compounds_most_severe_hazard_numerical_ratings = []
    compounds_hazard_sentences = []
    compounds_hazard_ratings = []
    compounds_hazard_colours = []
    compounds_exposure_potentials = []
    compounds_risk_ratings = []
    compounds_risk_colours = []
    # iterate through the compounds
    for compound_hazard_codes, physical_form in zip(
        compounds_hazard_codes_list, physical_forms_list
    ):
        # get all the data for the compound or skip if we don't have hazard codes and physical form data to use
        if not compound_hazard_codes or not physical_form:
            continue
        compound_data = get_single_compound_data(compound_hazard_codes, physical_form)
        # append the data to the relevant lists
        compounds_most_severe_hazard_numerical_ratings.append(
            compound_data["most_severe_hazard_numerical_rating"]
        )
        compounds_hazard_sentences.append(compound_data["hazard_sentence"])
        compounds_hazard_ratings.append(compound_data["hazard_ratings"])
        compounds_hazard_colours.append(compound_data["hazard_colour_code"])
        compounds_exposure_potentials.append(compound_data["exposure_potential"])
        compounds_risk_ratings.append(compound_data["risk_ratings"])
        compounds_risk_colours.append(compound_data["risk_colour_code"])

    return (
        compounds_most_severe_hazard_numerical_ratings,
        compounds_hazard_sentences,
        compounds_hazard_ratings,
        compounds_hazard_colours,
        compounds_exposure_potentials,
        compounds_risk_ratings,
        compounds_risk_colours,
    )


def get_single_compound_data(compound_hazard_codes: str, physical_form: str) -> Dict:
    """
    Gets the risk data for a single compound. To do so we must get hazard data to categorise the level of hazard
    The hazard category is combined with the exposure (determined by physical form) to give a risk score.

    Args:
        compound_hazard_codes - '-' delimited H codes. e.g., 'H301-H315' using the GHS hazard system
        physical_form - the physical form of the compound from the dropdown which relates to its exposure potential

    Returns:
        a dictionary containing:
            hazard_sentence - str: consisting of H-code and corresponding description
            hazard_ratings -list of strings ['H', 'L',...]
            most_serve_hazard_numerical_rating - int: A higher number means a more severe hazard
            hazard_colour_code - str: e.g., 'hazard-hazardous' indicates colouring of table cells
            exposure_potential - str: e.g., 'M'
            hazard_exposure_combination - str: e.g., 'VHM'
            risk - str: e.g., VH
            risk_colour_code - str: e.g.,'hazard-hazardous' indicates colouring of the table cells.

    """
    # first get the hazard data
    hazard_sentence, hazard_ratings = get_single_compound_hazard(compound_hazard_codes)
    # find the most severe hazard rating as int and str
    most_severe_hazard_numerical_rating = max(
        [numerical_rating_from_str[rating] for rating in hazard_ratings]
    )
    most_severe_hazard_rating = str_rating_from_number[
        most_severe_hazard_numerical_rating
    ]
    # Makes any very high (VH) hazards summary table cell red
    hazard_colour_code = (
        "hazard-hazardous"
        if most_severe_hazard_rating == "VH"
        else "hazard-reset-hazard"
    )
    # exposure data from the physical form. (volatile liquids have more exposure than dense solids)
    exposure_potential = physical_form_exposure_potential.get(physical_form)
    # hazard and exposure potential data are combined to make a risk rating for a compound
    hazard_exposure_combination = most_severe_hazard_rating + exposure_potential
    risk_rating = hazard_exposure_risk.get(hazard_exposure_combination)
    risk_colour_code = (
        "hazard-hazardous" if risk_rating == "VH" else "hazard-reset-hazard"
    )
    return {
        "hazard_sentence": hazard_sentence,
        "hazard_ratings": most_severe_hazard_rating,
        "most_severe_hazard_numerical_rating": most_severe_hazard_numerical_rating,
        "hazard_colour_code": hazard_colour_code,
        "exposure_potential": exposure_potential,
        "hazard_exposure_combination": hazard_exposure_combination,
        "risk_ratings": risk_rating,
        "risk_colour_code": risk_colour_code,
    }


def get_single_compound_hazard(compound_hazard_codes: str) -> Tuple[str, List[str]]:
    """
    For a single compound use the h codes, if present, to find and return the hazard sentences, and categories.

    Args:
        compound_hazard_codes - either delimited h codes or a str describing the lack of hazard data for this compound

    Returns:
        compound_hazard_sentence - combines the H-code and associated phrase
        compound_hazard_ratings - list of ['VH', 'L',...] low to very high hazard level of hazard codes for a compound
    """

    if (
        compound_hazard_codes == "Not Hazardous"  # solvent table value from chem21
        or compound_hazard_codes
        == "No hazard codes found"  # compound table value from pubchem
    ):
        # if there is no hazard data for the compound return these null values
        return compound_hazard_codes, ["L"]
    # if a string of hazard codes is present we first convert this to a list
    hazard_code_list = string_to_list(compound_hazard_codes)
    # we iterate through this list and for each hazard code we append the data to the appropriate list
    compound_hazard_sentence = ""
    compound_hazard_ratings = []
    for hazard_code in hazard_code_list:
        phrase, rating = get_data(hazard_code)
        rating = "M" if rating == "L" or rating == "" else rating
        compound_hazard_sentence += f"{hazard_code} {phrase}, "
        compound_hazard_ratings.append(rating)
    # post processing to remove delimiter from sentence string
    compound_hazard_sentence = compound_hazard_sentence[:-2]
    # we want the maximum risk rate from the compound_hazard_categories
    return compound_hazard_sentence, compound_hazard_ratings


def replace_substrings(input_string: str, replace_dict: dict) -> str:
    """
    Replace substrings in the input string using a dictionary of replacements.

    Args:
        input_string (str): The input string in which replacements should be made.
        replace_dict (dict): A dictionary where keys represent the substrings to be replaced
            and values represent the replacement strings.

    Returns:
        str: The input string with substrings replaced according to the dictionary.
    """
    for key, value in replace_dict.items():
        input_string = input_string.replace(key, value)
    return input_string


def get_data(h_code: str) -> Tuple[str, str]:
    """
    Gets the hazard data for a specific h_code e.g., H301 or H201

    Args:
        h_code

    Returns:
        the H-codes associated phrase
        the string rating of the hazard level
    """
    hazard_object = get(h_code)
    return (
        hazard_object.phrase,
        hazard_object.category if hazard_object.category is not None else "L",
    )


def string_to_list(compound_hazard_codes: str) -> List[str]:
    """
    Converts a string with variable delimiters into a list.

    Args:
        compound_hazard_codes - typically just '-' delimiter but sometimes others. sometimes in pubchem
        h codes are combined e.g., H301+H311 or other delimiters are used

    Returns:
        a list of hazard codes with each h code as an item in the list

    """
    compound_hazard_codes = replace_substrings(
        compound_hazard_codes, {"   ": ",", "-": ",", " + ": ",", "+": ","}
    )
    return sorted(list(set(compound_hazard_codes.split(","))))


def get(h_code: str) -> models.HazardCode:
    """
    Finds the database entry for a specific H-code

    Args:
        h_code - GHS H-code.
    Returns:
        a hazard code entry from the database containing additional hazard data.
    """
    return (
        db.session.query(models.HazardCode)
        .filter(models.HazardCode.code == h_code)
        .first()
    )
