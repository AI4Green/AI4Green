import json
from typing import Dict, List

from sources import auxiliary, models, services

empty_summary_table = json.dumps(
    {
        "summary_table_generated": False,
        "real_product_mass": "",
        "unreacted_reactant_mass": "",
        "polymer_mn": "",
        "polymer_mw": "",
        "polymer_dispersity": "",
        "polymer_mass_method": "-select-",
        "polymer_mass_calibration": "",
        "polymer_tg": "",
        "polymer_tm": "",
        "polymer_tc": "",
        "polymer_thermal_method": "-select-",
        "polymer_thermal_calibration": "",
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
        "reactant_primary_keys": [],
        "reagent_primary_keys": [],
        "solvent_primary_keys": [],
        "product_primary_keys": [],
    }
)


def get_request_data_from_keys(request_data: Dict, keys: List[str]) -> Dict:
    component_dict = {}
    for key in keys:
        snake_case_key = services.utils.camelCase_to_snake_case(key)
        component_dict[snake_case_key] = auxiliary.get_data(key, request_data)
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
        "reactantMolecularWeights",
        "reactantDensities",
        "reactantConcentrations",
        "reactantMns",
        "reactantEquivalents",
        "reactantAmounts",
        "roundedReactantAmounts",
        "reactantVolumes",
        "roundedReactantVolumes",
        "reactantMasses",
        "roundedReactantMasses",
        "reactantHazards",
        "reactantPhysicalForms",
        "reactantPrimaryKeys",
    ]
    reactant_dict = get_request_data_from_keys(request_data, keys)
    services.hazard_code.get_multiple_compounds_data(reactant_dict, "reactant")

    # Get sum values
    reactant_dict["reactant_mass_sum"] = float(request_data["reactantMassSum"])
    reactant_dict["reactant_molecular_weight_sum"] = float(
        request_data["reactantMolecularWeightSum"]
    )
    # Joining primary keys and getting smiles
    reactant_primary_keys_ls = auxiliary.get_data("reactantPrimaryKeys", request_data)
    reactant_dict["reactant_primary_keys_str"] = ", ".join(reactant_primary_keys_ls)
    polymer_indices = request_data["polymerIndices"]
    polymer_indices = json.loads(polymer_indices)
    reactant_dict["reactant_smiles_ls"] = services.all_compounds.get_smiles_list(
        reactant_primary_keys_ls, polymer_indices
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
        "reagentTableNumbers",
        "reagents",
        "reagentMolecularWeights",
        "reagentDensities",
        "reagentConcentrations",
        "reagentEquivalents",
        "reagentAmounts",
        "roundedReagentAmounts",
        "reagentVolumes",
        "roundedReagentVolumes",
        "reagentMasses",
        "roundedReagentMasses",
        "reagentHazards",
        "reagentPhysicalForms",
        "reagentPrimaryKeys",
    ]
    reagent_dict = get_request_data_from_keys(request_data, keys)
    services.hazard_code.get_multiple_compounds_data(reagent_dict, "reagent")

    # Get sum data
    reagent_dict["reagent_molecular_weight_sum"] = float(
        request_data["reagentMolecularWeightSum"]
    )
    # Getting primary keys and format
    reagent_primary_keys_ls = auxiliary.get_data("reagentPrimaryKeys", request_data)
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
        "solventTableNumbers",
        "solvents",
        "solventVolumes",
        "solventHazards",
        "solventPhysicalForms",
        "solventPrimaryKeys",
    ]
    solvent_dict = get_request_data_from_keys(request_data, solvent_data_keys)
    services.hazard_code.get_multiple_compounds_data(solvent_dict, "solvent")

    # Concurrency for if no solvents are chosen
    solvent_dict["number_of_solvents"] = request_data["numberOfSolvents"]
    if solvent_dict["number_of_solvents"] == "0":  # if no solvents have been chosen
        solvent_dict["number_of_solvents"] = 1  # then it shows only one empty cell
        solvent_dict["solvents"] = [" "]

    # handle solvent primary keys
    solvent_primary_keys_ls = auxiliary.get_data("solventPrimaryKeys", request_data)
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
        "productTableNumbers",
        "products",
        "productMasses",
        "roundedProductMasses",
        "productMolecularWeights",
        "productMns",
        "productEquivalents",
        "productHazards",
        "productPhysicalForms",
        "productPrimaryKeys",
    ]
    product_dict = get_request_data_from_keys(request_data, product_data_keys)
    services.hazard_code.get_multiple_compounds_data(product_dict, "product")

    # Primary keys and SMILES
    # Joining primary keys and getting smiles
    product_primary_keys_ls = auxiliary.get_data("productPrimaryKeys", request_data)

    polymer_indices = json.loads(request_data["polymerIndices"])
    polymer_indices = [
        (int(item) - int(request_data["numberOfReactants"])) for item in polymer_indices
    ]  # minus num reactants to update indexing

    product_dict["product_primary_keys_str"] = ", ".join(product_primary_keys_ls)
    product_dict["product_smiles_ls"] = services.all_compounds.get_smiles_list(
        product_primary_keys_ls, polymer_indices
    )
    # Table numbers
    product_table_numbers = auxiliary.get_data("productTableNumbers", request_data)
    product_dict["product_table_numbers"] = list(filter(None, product_table_numbers))
    product_dict["main_product_table_number"] = int(
        request_data["mainProductTableNumber"]
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
        "amount_unit": str(request_data["amountUnit"]),
        "volume_unit": str(request_data["volumeUnit"]),
        "mass_unit": str(request_data["massUnit"]),
        "solvent_volume_unit": str(request_data["solventVolumeUnit"]),
        "product_mass_unit": str(request_data["productMassUnit"]),
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
    # get a list of each hazard code in the reaction.
    risk_data["hazard_codes"] = list(
        set(
            reactants["reactant_hazard_codes"]
            + reagents["reagent_hazard_codes"]
            + solvents["solvent_hazard_codes"]
            + products["product_hazard_codes"]
        )
    )
    risk_data["toxicity_types"] = services.hazard_code.get_toxicities(
        risk_data["hazard_codes"]
    )
    return risk_data


def check_if_updated(
    reaction: models.Reaction,
    reactant_data: Dict,
    reagent_data: Dict,
    solvent_data: Dict,
    product_data: Dict,
) -> bool:
    """
    Checks if the reactants, reagents, solvents, or products in a reaction have changed.

    Args:
        reaction (models.Reaction): The reaction object containing summary table data
        reactant_data (Dict): Dictionary containing reactant information with 'reactant_primary_keys'
        reagent_data (Dict): Dictionary containing reagent information with 'reagent_primary_keys'
        solvent_data (Dict): Dictionary containing solvent information with 'solvent_primary_keys'
        product_data (Dict): Dictionary containing product information with 'product_primary_keys'

    Returns:
        bool: True if any component (reactants, reagents, solvents, or products) has changed,
              False if no changes detected

    Note:
        The function compares current primary keys with those stored in reaction.summary_table_data
    """
    summary_data = json.loads(reaction.summary_table_data)
    if reactant_data["reactant_primary_keys"] != summary_data.get(
        "reactant_primary_keys"
    ):
        return True
    if reagent_data["reagent_primary_keys"] != summary_data.get("reagent_primary_keys"):
        return True
    if solvent_data["solvent_primary_keys"] != summary_data.get("solvent_primary_keys"):
        return True
    if product_data["product_primary_keys"] != summary_data.get("product_primary_keys"):
        return True
    return False


def update_summary_keys(
    reaction: models.Reaction,
    reactant_data: Dict,
    reagent_data: Dict,
    solvent_data: Dict,
    product_data: Dict,
) -> None:
    """
    Updates the primary keys in a reaction's summary table data while preserving other fields.

    Args:
        reaction (models.Reaction): The reaction object to update with new summary data
        reactant_data (Dict): Dictionary containing reactant information with 'reactant_primary_keys'
        reagent_data (Dict): Dictionary containing reagent information with 'reagent_primary_keys'
        solvent_data (Dict): Dictionary containing solvent information with 'solvent_primary_keys'
        product_data (Dict): Dictionary containing product information with 'product_primary_keys'

    Returns:
        None: Updates the reaction.summary_table_data in place

    Note:
        The function modifies only the primary key fields in the existing summary_table_data,
        preserving all other existing keys and values
    """
    summary_data = json.loads(reaction.summary_table_data)

    summary_data["reactant_primary_keys"] = reactant_data["reactant_primary_keys"]
    summary_data["reagent_primary_keys"] = reagent_data["reagent_primary_keys"]
    summary_data["solvent_primary_keys"] = solvent_data["solvent_primary_keys"]
    summary_data["product_primary_keys"] = product_data["product_primary_keys"]

    reaction.summary_table_data = json.dumps(summary_data)
    reaction.update()
