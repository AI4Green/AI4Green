import re
from typing import Dict, List, Tuple

from flask import abort, request
from sources import auxiliary


def get_reactant_data() -> Dict:
    """
    Get reactant data from the request
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
    # Create a dictionary to store the data
    reactant_dict = {}
    # Iterate through each key and retrieve the data
    for key in keys:
        # Convert camelCase to snake_case
        snake_case_key = "".join(
            ["_" + c.lower() if c.isupper() else c for c in key]
        ).lstrip("_")
        reactant_dict[snake_case_key] = auxiliary.get_data(key)

    # Get sum values
    reactant_dict["reactant_mass_sum"] = float(request.form["reactantMassSum"])
    reactant_dict["reactant_molecular_weight_sum"] = float(
        request.form["reactantMolecularWeightSum"]
    )
    # Joining primary keys
    reactant_dict["reactant_primary_keys"] = ", ".join(reactant_dict["reactants"])
    return reactant_dict


def get_reagent_data() -> Dict:
    """
    Get reagent data from the request
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
    ]

    # Create a dictionary to store the reagent data
    reagent_dict = {}
    # Iterate through each key and retrieve the data
    for key in keys:
        # Convert camelCase to snake_case
        snake_case_key = "".join(
            ["_" + c.lower() if c.isupper() else c for c in key]
        ).lstrip("_")
        reagent_dict[snake_case_key] = auxiliary.get_data(key)

    # Get sum data
    reagent_dict["reagent_molecular_weight_sum"] = float(
        request.form["reagentMolecularWeightSum"]
    )
    # Getting primary keys and format
    reagent_primary_keys_ls = auxiliary.get_data("reagentPrimaryKeys")
    reagent_dict["reagent_primary_keys_str"] = ", ".join(reagent_primary_keys_ls)
    reagent_primary_keys_ls = [
        int(x) if x.isdigit() else reform_novel_compound_primary_key(x)
        for x in reagent_primary_keys_ls
        if x
    ]
    reagent_dict["reagent_primary_keys_list"] = reagent_primary_keys_ls
    return reagent_dict


def get_solvent_data() -> Dict:
    """
    Get solvent data from the request
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
    solvent_dict = {}

    # Iterate through each key and retrieve the data
    for key in solvent_data_keys:
        # Convert camelCase to snake_case
        snake_case_key = "".join(
            ["_" + c.lower() if c.isupper() else c for c in key]
        ).lstrip("_")
        solvent_dict[snake_case_key] = auxiliary.get_data(key)

    # Additional data
    solvent_dict["number_of_solvents"] = request.form["numberOfSolvents"]
    solvent_primary_keys_ls = auxiliary.get_data("solventPrimaryKeys")
    solvent_dict["solvent_primary_keys_str"] = ", ".join(solvent_primary_keys_ls)
    solvent_primary_keys_ls = [
        int(x) if x.isdigit() else reform_novel_compound_primary_key(x)
        for x in solvent_primary_keys_ls
        if x
    ]
    solvent_dict["solvent_primary_keys_list"] = solvent_primary_keys_ls
    return solvent_dict


# TODO
def make_rxn_file():
    return {}


def get_sustainability_metrics():
    return {}


def get_product_data() -> Dict:
    """
    Get product data from the request
    Returns:
         Dictionary with product data stored as a list with one entry per product or a single value.
    """
    product_data_keys = [
        "productTableNumbers",
        "products",
        "productMasses",
        "roundedProductMasses",
        "productMolecularWeights",
        "productHazards",
        "productPhysicalForms",
        "productPrimaryKeys",
    ]
    # Create a dictionary to store the product data
    product_dict = {}

    # Iterate through each key and retrieve the data
    for key in product_data_keys:
        # Convert camelCase to snake_case
        snake_case_key = "".join(
            ["_" + c.lower() if c.isupper() else c for c in key]
        ).lstrip("_")
        product_dict[snake_case_key] = auxiliary.get_data(key)

    # Primary keys
    product_primary_keys = auxiliary.get_data("productPrimaryKeys")
    product_dict["product_primary_keys"] = ", ".join(product_primary_keys)
    # Table numbers
    product_table_numbers = auxiliary.get_data("productTableNumbers")
    product_dict["product_table_numbers"] = list(filter(None, product_table_numbers))
    product_dict["main_product_table_number"] = int(
        request.form["mainProductTableNumber"]
    )

    for idx, product_table_num in enumerate(product_dict["product_table_numbers"]):
        if product_table_num == product_dict["main_product_table_number"]:
            product_dict["main_product_index"] = idx
            break
    return product_dict


def get_unit_data() -> Dict:
    """
    Gets a dictionary for the unit data from the reaction table.
    Returns:
        dictionary with unit data
    """
    # Gets selected units from the reaction table
    unit_dict = {
        "amount_unit": str(request.form["amountUnit"]),
        "volume_unit": str(request.form["volumeUnit"]),
        "mass_unit": str(request.form["massUnit"]),
        "solvent_volume_unit": str(request.form["solventVolumeUnit"]),
        "product_mass_unit": str(request.form["productMassUnit"]),
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
        reactant_data["reactantEquivalents"], reactant_data["reactantMasses"]
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
    if reagent_data.get("reagentTableNumbers"):
        for equivalents, hazard, name in zip(
            reagent_data["reagentEquivalents"],
            reagent_data["reagentHazards"],
            reagent_data["reagents"],
        ):
            if equivalents == "" or name == "" or hazard == "":
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
    if solvent_dict.get("solventTableNumbers"):
        for sol, vol in zip(solvent_dict["solvents"], solvent_dict["solventVolumes"]):
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
            reactant_data["reactantPhysicalForms"]
            + reagent_data.get("reagentPhysicalForms", [])
            + solvent_dict.get("solventPhysicalForms", [])
            + product_dict["productPhysicalForms"]
        ),
    ]

    for check_result in checks:
        if check_result != "check successful":
            return check_result
    return "checks successful"


def reform_novel_compound_primary_key(primary_key: str) -> Tuple:
    """
    Converts a novel primary key to a tuple from the string returned from the frontend HTML

    Args:
        primary_key - the primary key as a string. e.g., ('pixie dust', 1)

    Returns:
        A tuple of (compound_name, workbook_id)
    """
    if len(primary_key) > 350:
        abort(
            413
        )  # content too large. Exceeds max workbook name length + max novel compound name length
    compound_name = re.search(r"\('([^']*)', \d", primary_key).group(1)
    workbook_id = int(re.search(r"', (\d+)", primary_key).group(1))
    return compound_name, workbook_id
