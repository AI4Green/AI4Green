# from typing import List
#
# from sources import models
# from sources.extensions import db
#
# #
# #         compound_hazard_sentences,
# #         compound_hazard_ratings,
# #         compound_hazard_colors,
# #         compound_risk_colors,
# #         compound_exposure_potentials,
# #         compound_risk_ratings,
# #         total_rate,
#
#
# category_rate = {
#     "": 0,
#     "L": 1,
#     "M": 2,
#     "H": 3,
#     "VH": 4,
# }  # 'risk category - rate' dictionary
#
# category_rate_rev = {
#     "": 0,
#     "L": 1,
#     "M": 2,
#     "H": 3,
#     "VH": 4,
# }  # 'risk category - rate' dictionary
#
# physical_form_exposure_potential = {  # physical form - risk category' dictionary
#     "-select-": "M",
#     "Dense solid": "L",
#     "Non-volatile liquid": "L",
#     "Unknown": "M",
#     "Dusty Solid": "M",
#     "Lyophilised solid": "M",
#     "Volatile liquid": "M",
#     "Gas": "H",
#     "Highly volatile liquid": "H",
#     "Aerosol": "H",
#     "Solution that promotes skin absorption": "H",
# }
# hazard_exposure_risk = {  # 'combination of risks - risk category' dictionary
#     "VHL": "H",
#     "VHM": "VH",
#     "VHH": "VH",
#     "HL": "M",
#     "HM": "H",
#     "HH": "H",
#     "ML": "L",
#     "MM": "M",
#     "MH": "H",
#     "LL": "L",
#     "LM": "L",
#     "LH": "M",
# }
#
#
# def get_multiple_compound_hazard_data(
#     compounds_hazard_codes_list: List[str], physical_forms_list: List[str]
# ):
#     """
#     For a list of compounds (e.g., products) this function returns their hazard data
#
#     This function gets compounds, their hazards and physical forms, and returns
#     hazard sentences, hazard ratings, hazard and risk color codes, exposure potentials,
#     risk ratings, and total hazard rates
#
#
#     Args:
#         compounds_hazard_codes_list example: ['H201-H302', 'H401-H310']
#
#     """
#
#     # initiate the list variables to be returned
#     total_rate = []
#     compounds_hazard_sentences = []
#     compounds_hazard_ratings = []
#     compounds_hazard_colours = []
#     compounds_exposure_potentials = []
#     compounds_risk_ratings = []
#     compounds_risk_colours = []
#     # iterate through the compounds
#     for compound_hazard_codes, physical_form in zip(
#         compounds_hazard_codes_list, physical_forms_list
#     ):
#         (
#             highest_rated_hazard,
#             hazard_sentence,
#             hazard_rating,
#             hazard_colour,
#             exposure_potential,
#             risk_rating,
#             risk_colour,
#         ) = get_single_compound_risk_data(compound_hazard_codes, physical_form)
#         total_rate.append(highest_rated_hazard)
#         compounds_hazard_sentences.append(hazard_sentence)
#         compounds_hazard_ratings.append(hazard_rating)
#         compounds_hazard_colours.append(hazard_colour)
#         compounds_exposure_potentials.append(exposure_potential)
#         compounds_risk_ratings.append(risk_rating)
#         compounds_risk_colours.append(risk_colour)
#
#     return (
#         compounds_hazard_sentences,
#         compounds_hazard_ratings,
#         compounds_hazard_colours,
#         compounds_risk_colours,
#         compounds_exposure_potentials,
#         compounds_risk_ratings,
#         total_rate,
#     )
#
#
# def get_single_compound_risk_data(compound_hazard_codes, physical_form):
#     """
#     Gets the risk data for a single compound. To do so we must get hazard data to categorise the level of hazard
#     The hazard category is combined with the exposure (determined by physical form) to give a risk score.
#
#     :param compound_hazard_codes:
#     :param physical_form:
#     :return:
#     """
#     # first get the hazard data
#     hazard_sentence, highest_hazard_rating = get_single_compound_hazard(
#         compound_hazard_codes
#     )
#     compounds_hazard_sentences.append(hazard_sentence)
#     total_rate.append(highest_hazard_rating)
#     hazard_rating = category_rate_rev[highest_hazard_rating]
#     compounds_hazard_ratings.append(hazard_rating)
#     # in format 'hazard-...' as in models.User colours. Makes any very high (VH) hazards summary table cell red
#     hazard_colour_code = (
#         "hazard-hazardous" if hazard_rating == "VH" else "hazard-reset-hazard"
#     )
#     compounds_hazard_colours.append(hazard_colour_code)
#     # exposure and risk
#     exposure_potential = physical_form_exposure_potential.get(physical_form)
#     hazard_exposure_combination = hazard_rating + exposure_potential
#     risk_rating = hazard_exposure_risk.get(hazard_exposure_combination)
#     compounds_risk_ratings.append(risk_rating)
#     risk_colour_code = (
#         "hazard-hazardous" if risk_rating == "VH" else "hazard-reset-hazard"
#     )
#     compounds_risk_colours.append(risk_colour_code)
#
#
# def get_compound_exposure_potentials(component):
#     compound_exposure_potentials = []  # list of exposure potentials
#     for compound_physical_form in physical_forms_list:
#         # converts physical forms to exposure potentials
#         compound_exposure_potential = physical_form_exposure_potential.get(
#             compound_physical_form
#         )
#         compound_exposure_potentials.append(compound_exposure_potential)
#     compound_hazard_exposures = [
#         compound_hazard_rating + compound_exposure_potential
#         for (compound_hazard_rating, compound_exposure_potential) in zip(
#             compound_hazard_ratings, compound_exposure_potentials
#         )
#     ]
#     # combines hazard ratings and exposure potentials
#     compound_risk_ratings = []  # the risk rating list
#     compound_risk_colors = []  # the risk colour list
#     for compound_hazard_exposure in compound_hazard_exposures:
#         # converts hazard exposures to risk ratings and assigns color codes to them
#         compound_risk_rating = hazard_exposure_risk.get(compound_hazard_exposure)
#         compound_risk_ratings.append(compound_risk_rating)
#         compound_risk_color = (
#             "hazard-hazardous"
#             if compound_risk_rating == "VH"
#             else "hazard-reset-hazard"
#         )
#         compound_risk_colors.append(compound_risk_color)
#
#
# def replace_substrings(input_string: str, replace_dict: dict) -> str:
#     """
#     Replace substrings in the input string using a dictionary of replacements.
#
#     Args:
#         input_string (str): The input string in which replacements should be made.
#         replace_dict (dict): A dictionary where keys represent the substrings to be replaced
#             and values represent the replacement strings.
#
#     Returns:
#         str: The input string with substrings replaced according to the dictionary.
#     """
#     for key, value in replace_dict.items():
#         input_string = input_string.replace(key, value)
#     return input_string
#
#
# # def get_compounds_hazard_data(compound_hazard_codes):
# #     if compound_hazard_codes in {"Not Hazardous", "No hazard codes found"}:
# #         compound_sentence = 'compound_hazard'
# #         compound_hazard_rating = 'L'
# #     compound_hazard_codes_list = make_hazard_code_list(compound_hazard_codes)
# #     compound_hazard_sentence = ''
# #     for hazard_code in compound_hazard_codes_list:
# #         hazard_code_data = get_hazard_data_from_h_code(hazard_code)
# #         compound_hazard_sentence += f'{hazard_code} {hazard_code_data.phrase}, '
# #
# #
# #         # hazard_phrase = (
# #         #     hazard_code_object.phrase
# #         # )  # appends the risk rate to its list
# #         # compound_hazard_sentence += (
# #         #         hazard_code + " " + hazard_phrase + ", "
# #         # )  # the hazard sentence
# #         hazard_category = category_rate.get(hazard_code_object.category)
# #         if hazard_category is not None:
# #             compound_hazard_categories.append(hazard_category)
# #
# #     return 'x'
#
#
# def get_single_compound_hazard(compound_hazard_codes: str):
#     """
#     For a single compound interpret the hazard data present - either delimtied h codes or a string describing the lack
#     of hazard data for this compound.
#     Then use these identifiers to retrieve the detailed hazard phrases, and categories
#
#     :param compound_hazard_codes:
#     :return:
#     """
#
#     if (
#         compound_hazard_codes == "Not Hazardous"
#         or compound_hazard_codes == "No hazard codes found"
#     ):
#         return get_null_data_values(compound_hazard_codes)
#     # if a string of hazard codes is present we first convert this to a list
#     hazard_code_list = make_hazard_code_list(compound_hazard_codes)
#     # we iterate through this list and for each hazard code we append the data to the appropriate list
#     compound_hazard_sentence = ""
#     compound_hazard_categories = []
#     for hazard_code in hazard_code_list:
#         phrase, category = get_data_from_h_code(hazard_code)
#         compound_hazard_sentence += f"{hazard_code} {phrase}, "
#         # compound_rate = category_rate.get(compound.category)
#         # in the database they are described as L, M, H, VH. We change these to corresponding numbers
#         compound_hazard_categories.append(category_rate[category])
#     # post processing to remove delimiter from sentence string
#     compound_hazard_sentence = compound_hazard_sentence[:-2]
#     # we want the maximum risk rate from the compound_hazard_categories
#     max_compound_rate = int(max(compound_hazard_categories))
#     return compound_hazard_sentence, max_compound_rate
#
#
# def get_null_data_values(null_data_statement: str):
#     """
#     Gets the sentence and rating for compounds with no hazard data in the database
#
#     """
#     compound_hazard_sentence = null_data_statement
#     compound_hazard_rating = "L"
#     return compound_hazard_sentence, compound_hazard_rating
#
#
# def get_data_from_h_code(h_code: str):
#     """
#     Gets the hazard data for a specific h_code e.g., H301 or H201
#
#     Args
#         h_code
#
#     Returns
#     """
#     hazard_object = get_hazardcode_database_object(h_code)
#     return (
#         hazard_object.phrase,
#         hazard_object.category if hazard_object.category is not None else "L",
#     )
#
#
# def make_hazard_code_list(compound_hazard_codes: str):
#     """
#     Converts a string with variable delimiters into a list.
#     """
#     compound_hazard_codes = replace_substrings(
#         compound_hazard_codes, {"   ": ",", "-": ",", " + ": ",", "+": ","}
#     )
#     return sorted(list(set(compound_hazard_codes.split(","))))
#
#
# def get_hazardcode_database_object(h_code: str) -> models.HazardCode:
#     return (
#         db.session.query(models.HazardCode)
#         .filter(models.HazardCode.code == h_code)
#         .first()
#     )
