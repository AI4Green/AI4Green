from typing import Dict, List

from flask import render_template
from sources import auxiliary, services


def get_request_data_from_keys(request_data: Dict, keys: List[str]) -> Dict:
    component_dict = {}
    for key in keys:
        # snake_case_key = services.utils.camelCase_to_snake_case(key)
        component_dict[key] = auxiliary.get_data(key, request_data)
    return component_dict


def get_reactant_data(request_data: Dict) -> Dict:
    """
    Get reactant data from the request
    Args:
        request_data is the dictionary from the flask.request.form/json

    Returns:
         Dictionary with reactant data stored as a list with one entry per reactant and single floats for sum values.
    """
    keys = [
        "reactants",
        "reactant_molecular_weights",
        "reactant_densities",
        "reactant_concentrations",
        "reactant_equivalents",
        "reactant_amounts",
        "rounded_reactant_amounts",
        "reactant_volumes",
        "rounded_reactant_volumes",
        "reactant_masses",
        "rounded_reactant_masses",
        "reactant_hazards",
        "reactant_physical_forms",
    ]

    reactant_dict = get_request_data_from_keys(request_data, keys)
    services.hazard_code.get_multiple_compounds_data(reactant_dict, "reactant")

    # Get sum values
    reactant_dict["reactant_mass_sum"] = float(request_data["reactant_mass_sum"])
    reactant_dict["reactant_molecular_weight_sum"] = float(
        request_data["reactant_molecular_weight_sum"]
    )
    # Joining primary keys and getting smiles
    reactant_primary_keys_ls = auxiliary.get_data("reactant_primary_keys", request_data)
    reactant_dict["reactant_primary_keys_str"] = ", ".join(reactant_primary_keys_ls)
    reactant_dict["reactant_smiles_ls"] = services.all_compounds.get_smiles_list(
        reactant_primary_keys_ls
    )
    return reactant_dict


def get_reagent_data(request_data: Dict) -> Dict:
    """
    Get reagent data from the request
    Args:
        request_data is the dictionary from the flask.request.form/json

    Returns:
         Dictionary with reagent data stored as a list with one entry per reagent and single floats for sum values.
    """
    keys = [
        "reagent_table_numbers",
        "reagents",
        "reagent_molecular_weights",
        "reagent_densities",
        "reagent_concentrations",
        "reagent_equivalents",
        "reagent_amounts",
        "rounded_reagent_amounts",
        "reagent_volumes",
        "rounded_reagent_volumes",
        "reagent_masses",
        "rounded_reagent_masses",
        "reagent_hazards",
        "reagent_physical_forms",
    ]

    reagent_dict = get_request_data_from_keys(request_data, keys)
    services.hazard_code.get_multiple_compounds_data(reagent_dict, "reagent")

    # Get sum data
    reagent_dict["reagent_molecular_weight_sum"] = float(
        request_data["reagent_molecular_weight_sum"]
    )
    # Getting primary keys and format
    reagent_primary_keys_ls = auxiliary.get_data("reagent_primary_keys", request_data)
    reagent_dict["reagent_primary_keys_str"] = ", ".join(reagent_primary_keys_ls)
    reagent_dict["reagent_primary_keys_list"] = reagent_primary_keys_ls
    reagent_dict["reagent_smiles_ls"] = services.all_compounds.get_smiles_list(
        reagent_primary_keys_ls
    )
    return reagent_dict


def get_solvent_data(request_data: Dict) -> Dict:
    """
    Get solvent data from the request
    Args:
        request_data is the dictionary from the flask.request.form/json

    Returns:
         Dictionary with solvent data stored as a list with one entry per solvent or a single value.
    """
    solvent_data_keys = [
        "solvent_table_numbers",
        "solvents",
        "solvent_volumes",
        "solvent_hazards",
        "solvent_physical_forms",
    ]
    solvent_dict = get_request_data_from_keys(request_data, solvent_data_keys)
    services.hazard_code.get_multiple_compounds_data(solvent_dict, "solvent")

    # Concurrency for if no solvents are chosen
    solvent_dict["number_of_solvents"] = request_data["number_of_solvents"]
    if solvent_dict["number_of_solvents"] == "0":  # if no solvents have been chosen
        solvent_dict["number_of_solvents"] = 1  # then it shows only one empty cell
        solvent_dict["solvents"] = [" "]

    # handle solvent primary keys
    solvent_primary_keys_ls = auxiliary.get_data("solvent_primary_keys", request_data)
    solvent_dict["solvent_primary_keys_str"] = ", ".join(solvent_primary_keys_ls)
    solvent_dict["solvent_primary_keys_list"] = solvent_primary_keys_ls
    # get the smiles from the primary keys
    solvent_dict["solvent_smiles_ls"] = services.all_compounds.get_smiles_list(
        solvent_primary_keys_ls
    )
    return solvent_dict


def get_product_data(request_data: Dict) -> Dict:
    """
    Get product data from the request
    Args:
        request_data is the dictionary from the flask.request.form/json

    Returns:
         Dictionary with product data stored as a list with one entry per product or a single value.
    """
    product_data_keys = [
        "product_table_numbers",
        "products",
        "product_masses",
        "rounded_product_masses",
        "product_molecular_weights",
        "product_hazards",
        "product_physical_forms",
    ]
    product_dict = get_request_data_from_keys(request_data, product_data_keys)
    services.hazard_code.get_multiple_compounds_data(product_dict, "product")

    # Primary keys and SMILES
    # Joining primary keys and getting smiles
    product_primary_keys_ls = auxiliary.get_data("product_primary_keys", request_data)
    product_dict["product_primary_keys_str"] = ", ".join(product_primary_keys_ls)
    product_dict["product_smiles_ls"] = services.all_compounds.get_smiles_list(
        product_primary_keys_ls
    )
    # Table numbers
    product_table_numbers = auxiliary.get_data("product_table_numbers", request_data)
    product_dict["product_table_numbers"] = list(filter(None, product_table_numbers))
    product_dict["main_product_table_number"] = int(
        request_data["main_product_table_number"]
    )

    for idx, product_table_num in enumerate(product_dict["product_table_numbers"]):
        if int(product_table_num) == product_dict["main_product_table_number"]:
            product_dict["main_product_index"] = idx
            break
    return product_dict


def get_unit_data(request_data: Dict) -> Dict:
    """
    Gets a dictionary for the unit data from the reaction table.
    Args:
        request_data is the dictionary from the flask.request.form/json

    Returns:
        dictionary with unit data
    """
    # Gets selected units from the reaction table
    unit_dict = {
        "amount_unit": str(request_data["amount_unit"]),
        "volume_unit": str(request_data["volume_unit"]),
        "mass_unit": str(request_data["mass_unit"]),
        "solvent_volume_unit": str(request_data["solvent_volume_unit"]),
        "product_mass_unit": str(request_data["product_mass_unit"]),
    }
    return unit_dict


def check_reactant_data(reactant_data: Dict) -> str:
    """
    Check if all necessary reactant data is present.

    Args:
        reactant_data (Dict): Dictionary containing reactant data.

    Returns:
        str: Returns a message indicating missing data or "check successful" if all data is present.
    """
    for equivalents, mass in zip(
        reactant_data["reactant_equivalents"], reactant_data["reactant_masses"]
    ):
        if mass == "" or mass == 0 or equivalents == "":
            return "Ensure you have entered all the necessary reactant information!"
    return "check successful"


def check_reagent_data(reagent_data: Dict) -> str:
    """
    Check if all necessary reagent data is present.

    Args:
        reagent_data (Dict): Dictionary containing reagent data.

    Returns:
        str: Returns a message indicating missing data or "check successful" if all data is present.
    """
    if reagent_data.get("reagents")[0]:
        for equivalents, name in zip(
            reagent_data["reagent_equivalents"],
            reagent_data["reagents"],
        ):
            if equivalents == "" or name == "":
                return "Ensure you have entered all the necessary reagent information!"
    return "check successful"


def check_solvent_data(solvent_dict: Dict) -> str:
    """
    Check if all necessary solvent data is present.

    Args:
        solvent_dict (Dict): Dictionary containing solvent data.

    Returns:
        str: Returns a message indicating missing data or "check successful" if all data is present.
    """
    if solvent_dict.get("solvents")[0] and solvent_dict.get("solvents") != [" "]:
        for sol, vol in zip(solvent_dict["solvents"], solvent_dict["solvent_volumes"]):
            if vol == "" or vol == 0 or sol == "" or sol == 0:
                return "Ensure you have entered all the necessary solvent information!"
    return "check successful"


def check_physical_forms(components: List[str]) -> str:
    """
    Check if all necessary physical form data is present.

    Args:
        components (List[str]): List containing physical form data.

    Returns:
        str: Returns a message indicating missing data or "check successful" if all data is present.
    """
    for phys_form in components:
        if phys_form == "-select-":
            return (
                "Ensure you have entered all the necessary physical form information!"
            )
    return "check successful"


def check_required_data_is_present(
    reactant_data: Dict, reagent_data: Dict, solvent_dict: Dict, product_dict: Dict
) -> str:
    """
    Check if all required data that the user must enter is present.

    Args:
        reactant_data (Dict): Dictionary containing reactant data.
        reagent_data (Dict): Dictionary containing reagent data.
        solvent_dict (Dict): Dictionary containing solvent data.
        product_dict (Dict): Dictionary containing product data.

    Returns:
        str: Returns a message indicating missing data or "checks successful" if all data is present.
    """
    # user does not enter data about product at this stage other than physical form.
    checks = [
        check_reactant_data(reactant_data),
        check_reagent_data(reagent_data),
        check_solvent_data(solvent_dict),
        check_physical_forms(
            reactant_data["reactant_physical_forms"]
            + reagent_data.get("reagent_physical_forms", [])
            + solvent_dict.get("solvent_physical_forms", [])
            + product_dict["product_physical_forms"]
        ),
    ]

    for check_result in checks:
        if check_result != "check successful":
            return check_result
    return "checks successful"


def get_risk_data(
    reactants: Dict, reagents: Dict, solvents: Dict, products: Dict
) -> Dict:
    category_rate = {
        "": 0,
        "L": 1,
        "M": 2,
        "H": 3,
        "VH": 4,
    }

    risk_data = {}
    most_severe_hazard_numerical_rating = (
        reactants["reactant_most_severe_hazard_numerical_ratings"]
        + reagents["reagent_most_severe_hazard_numerical_ratings"]
        + solvents["solvent_most_severe_hazard_numerical_ratings"]
        + products["product_most_severe_hazard_numerical_ratings"]
    )
    max_most_severe_hazard_numerical_rating = int(
        max(most_severe_hazard_numerical_rating)
    )  # max total hazard rate
    risk_data["risk_rating"] = list(category_rate.keys())[
        max_most_severe_hazard_numerical_rating
    ]  # resulting hazard rating
    risk_data["risk_colour"] = (
        "hazard-hazardous"
        if risk_data["risk_rating"] == "VH"
        else "hazard-reset-hazard"
    )  # colour code for the hazard rating
    return risk_data


def render_summary_template(summary_data: Dict):
    return render_template(
        "_summary_table.html",
        amount_unit=summary_data["amount_unit"],
        volume_unit=summary_data["volume_unit"],
        mass_unit=summary_data["mass_unit"],
        solvent_volume_unit=summary_data["solvent_volume_unit"],
        product_mass_unit=summary_data["product_mass_unit"],
        reactants=summary_data["reactants"],
        reactant_primary_keys=summary_data["reactant_primary_keys_str"],
        reagent_primary_keys=summary_data["reagent_primary_keys_str"],
        reagents=summary_data["reagents"],
        reagent_table_numbers=summary_data["reagent_table_numbers"],
        reagent_molecular_weights=summary_data["reagent_molecular_weights"],
        reagent_densities=summary_data["reagent_densities"],
        reagent_concentrations=summary_data["reagent_concentrations"],
        reagent_equivalents=summary_data["reagent_equivalents"],
        reagent_hazards=summary_data["reagent_hazards"],
        reagent_amounts=summary_data["reagent_amounts"],
        rounded_reagent_amounts=summary_data["rounded_reagent_amounts"],
        reagent_volumes=summary_data["reagent_volumes"],
        rounded_reagent_volumes=summary_data["rounded_reagent_volumes"],
        reagent_masses=summary_data["reagent_masses"],
        rounded_reagent_masses=summary_data["rounded_reagent_masses"],
        solvents=summary_data["solvents"],
        solvent_volumes=summary_data["solvent_volumes"],
        solvent_table_numbers=summary_data["solvent_table_numbers"],
        solvent_flags=summary_data["solvent_flags"],
        products=summary_data["products"],
        product_table_numbers=summary_data["product_table_numbers"],
        reactant_molecular_weights=summary_data["reactant_molecular_weights"],
        reactant_densities=summary_data["reactant_densities"],
        reactant_concentrations=summary_data["reactant_concentrations"],
        reactant_equivalents=summary_data["reactant_equivalents"],
        reactant_amounts=summary_data["reactant_amounts"],
        rounded_reactant_amounts=summary_data["rounded_reactant_amounts"],
        reactant_volumes=summary_data["reactant_volumes"],
        rounded_reactant_volumes=summary_data["rounded_reactant_volumes"],
        reactant_masses=summary_data["reactant_masses"],
        rounded_reactant_masses=summary_data["rounded_reactant_masses"],
        product_primary_keys=summary_data["product_primary_keys_str"],
        main_product_table_number=summary_data["main_product_table_number"],
        main_product_index=summary_data["main_product_index"],
        product_molecular_weights=summary_data["product_molecular_weights"],
        product_masses=summary_data["product_masses"],
        rounded_product_masses=summary_data["rounded_product_masses"],
        ae=summary_data["ae"],
        ae_flag=summary_data["ae_flag"],
        element_sustainability=summary_data["element_sustainability"],
        element_sustainability_flag=summary_data["element_sustainability_flag"],
        reactant_hazard_sentences=summary_data["reactant_hazard_sentences"],
        reactant_hazard_ratings=summary_data["reactant_hazard_ratings"],
        reactant_hazard_colors=summary_data["reactant_hazard_colours"],
        reactant_risk_colors=summary_data["reactant_risk_colours"],
        reactant_exposure_potentials=summary_data["reactant_exposure_potentials"],
        reactant_risk_ratings=summary_data["reactant_risk_ratings"],
        reagent_hazard_sentences=summary_data["reagent_hazard_sentences"],
        reagent_hazard_ratings=summary_data["reagent_hazard_ratings"],
        reagent_hazard_colors=summary_data["reagent_hazard_colours"],
        reagent_risk_colors=summary_data["reagent_risk_colours"],
        reagent_exposure_potentials=summary_data["reagent_exposure_potentials"],
        reagent_risk_ratings=summary_data["reagent_risk_ratings"],
        solvent_primary_keys=summary_data["solvent_primary_keys_str"],
        solvent_hazard_sentences=summary_data["solvent_hazard_sentences"],
        solvent_hazard_ratings=summary_data["solvent_hazard_ratings"],
        solvent_exposure_potentials=summary_data["solvent_exposure_potentials"],
        solvent_risk_ratings=summary_data["solvent_risk_ratings"],
        solvent_hazard_colors=summary_data["solvent_hazard_colours"],
        solvent_risk_colors=summary_data["solvent_risk_colours"],
        product_hazard_sentences=summary_data["product_hazard_sentences"],
        product_hazard_ratings=summary_data["product_hazard_ratings"],
        product_exposure_potentials=summary_data["product_exposure_potentials"],
        product_risk_ratings=summary_data["product_risk_ratings"],
        product_hazard_colors=summary_data["product_hazard_colours"],
        product_risk_colors=summary_data["product_risk_colours"],
        risk_rating=summary_data["risk_rating"],
        risk_color=summary_data["risk_colour"],
        number_of_solvents=summary_data["number_of_solvents"],
    )
