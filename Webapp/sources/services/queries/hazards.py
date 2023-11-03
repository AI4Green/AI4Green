from typing import List

from sources import models
from sources.extensions import db

#
#         compound_hazard_sentences,
#         compound_hazard_ratings,
#         compound_hazard_colors,
#         compound_risk_colors,
#         compound_exposure_potentials,
#         compound_risk_ratings,
#         total_rate,


def get_compound_hazard_data(
    component_hazard_codes_list: List[str], component_physical_forms_list: List[str]
):
    """
    For a list of compounds (e.g., products) this function returns their hazard data

    This function gets compounds, their hazards and physical forms, and returns
    hazard sentences, hazard ratings, hazard and risk color codes, exposure potentials,
    risk ratings, and total hazard rates


    Args:
        component_hazard_codes_list example: ['H201-H302', 'H401-H310']

    TODO: Refactor this function.
    """

    (
        category_rate,
        physical_form_exposure_potential,
        hazard_exposure_risk,
    ) = get_hazard_dictionaries()
    #
    total_rate = [0]  # the list of maximum hazard rates
    compound_hazard_sentences = []  # the list of compound hazard sentences
    compound_hazard_ratings = []  # the list of compound hazard ratings
    compound_hazard_colors = []  # list of compound hazard colors

    for compound_hazard_codes in component_hazard_codes_list:
        compound_hazard_categories = [0]
        # sentences, ratings, colours = get_compounds_hazard_data(compound_hazard_codes)
        if compound_hazard_codes in {"Not Hazardous", "No hazard codes found"}:
            compound_hazard_sentence = compound_hazard_codes
            compound_hazard_rating = "L"
        compound_hazard_codes_list = make_hazard_code_list(compound_hazard_codes)
        compound_hazard_sentence = ""

        for hazard_code in compound_hazard_codes_list:
            hazard_code_object = get_hazard_data_from_h_code(hazard_code)
            compound_hazard_sentence += f"{hazard_code} {hazard_code_object.phrase}, "
            hazard_category = category_rate.get(hazard_code_object.category)
            if hazard_category is not None:
                compound_hazard_categories.append(hazard_category)

        # compound_hazard_sentence = ""  # joins hazard codes and phrases
        # for hazard_code in compound_hazard_codes:
        #     """We use a hazard code to query the hazard database
        #     and return a corresponding hazard phrase and risk category"""
        #     hazard_code_object = (
        #         db.session.query(models.HazardCode)
        #         .filter(models.HazardCode.code == hazard_code)
        #         .first()
        #     )
        #     print(hazard_code, hazard_code_object)
        #     hazard_phrase = (
        #         hazard_code_object.phrase
        #     )  # appends the risk rate to its list
        #     compound_hazard_sentence += (
        #         hazard_code + " " + hazard_phrase + ", "
        #     )  # the hazard sentence
        #     hazard_category = category_rate.get(hazard_code_object.category)

        compound_hazard_sentence = compound_hazard_sentence[
            :-2
        ]  # removes the comma and space at the end
        max_compound_rate = int(max(compound_hazard_categories))  # maximum risk rate
        total_rate.append(
            max_compound_rate
        )  # appends the maximum risk rate to their list

        compound_hazard_rating = (
            list(category_rate.keys())[max_compound_rate]
            if max_compound_rate > 0
            else "M"
        )  # converts the risk rates to hazard ratings
        # if the maximum risk rate is 0, then the hazard rating is M - medium
        # else:  # if the compound is not hazardous,
        #     compound_hazard_sentence = (
        #         compound_hazard  # then the hazard sentence is 'Not hazardous'
        #     )

        compound_hazard_sentences.append(
            compound_hazard_sentence
        )  # appends the hazard sentence to their list
        compound_hazard_ratings.append(
            compound_hazard_rating
        )  # appends the hazard rating to their list
        compound_hazard_color = (
            "hazard-hazardous"
            if compound_hazard_rating == "VH"
            else "hazard-reset-hazard"
        )  # assigns color codes
        compound_hazard_colors.append(
            compound_hazard_color
        )  # appending the table cell colour to their list
    compound_exposure_potentials = []  # list of exposure potentials
    for compound_physical_form in component_physical_forms_list:
        # converts physical forms to exposure potentials
        compound_exposure_potential = physical_form_exposure_potential.get(
            compound_physical_form
        )
        compound_exposure_potentials.append(compound_exposure_potential)
    compound_hazard_exposures = [
        compound_hazard_rating + compound_exposure_potential
        for (compound_hazard_rating, compound_exposure_potential) in zip(
            compound_hazard_ratings, compound_exposure_potentials
        )
    ]
    # combines hazard ratings and exposure potentials
    compound_risk_ratings = []  # the risk rating list
    compound_risk_colors = []  # the risk colour list
    for compound_hazard_exposure in compound_hazard_exposures:
        # converts hazard exposures to risk ratings and assigns color codes to them
        compound_risk_rating = hazard_exposure_risk.get(compound_hazard_exposure)
        compound_risk_ratings.append(compound_risk_rating)
        compound_risk_color = (
            "hazard-hazardous"
            if compound_risk_rating == "VH"
            else "hazard-reset-hazard"
        )
        compound_risk_colors.append(compound_risk_color)

    return (
        compound_hazard_sentences,
        compound_hazard_ratings,
        compound_hazard_colors,
        compound_risk_colors,
        compound_exposure_potentials,
        compound_risk_ratings,
        total_rate,
    )


def get_hazard_dictionaries():
    category_rate = {
        "": 0,
        "L": 1,
        "M": 2,
        "H": 3,
        "VH": 4,
    }  # 'risk category - rate' dictionary

    physical_form_exposure_potential = {  # physical form - risk category' dictionary
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
    hazard_exposure_risk = {  # 'combination of risks - risk category' dictionary
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
    return category_rate, physical_form_exposure_potential, hazard_exposure_risk


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


# def get_compounds_hazard_data(compound_hazard_codes):
#     if compound_hazard_codes in {"Not Hazardous", "No hazard codes found"}:
#         compound_sentence = 'compound_hazard'
#         compound_hazard_rating = 'L'
#     compound_hazard_codes_list = make_hazard_code_list(compound_hazard_codes)
#     compound_hazard_sentence = ''
#     for hazard_code in compound_hazard_codes_list:
#         hazard_code_data = get_hazard_data_from_h_code(hazard_code)
#         compound_hazard_sentence += f'{hazard_code} {hazard_code_data.phrase}, '
#
#
#         # hazard_phrase = (
#         #     hazard_code_object.phrase
#         # )  # appends the risk rate to its list
#         # compound_hazard_sentence += (
#         #         hazard_code + " " + hazard_phrase + ", "
#         # )  # the hazard sentence
#         hazard_category = category_rate.get(hazard_code_object.category)
#         if hazard_category is not None:
#             compound_hazard_categories.append(hazard_category)
#
#     return 'x'


def make_hazard_code_list(compound_hazard_codes: str):
    """
    Converts a string with variable delimiters into a list.
    """
    compound_hazard_codes = replace_substrings(
        compound_hazard_codes, {"   ": ",", "-": ",", " + ": ",", "+": ","}
    )
    return sorted(list(set(compound_hazard_codes.split(","))))


def get_hazard_data_from_h_code(h_code: str) -> models.HazardCode:
    return (
        db.session.query(models.HazardCode)
        .filter(models.HazardCode.code == h_code)
        .first()
    )
