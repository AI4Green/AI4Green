import re
from typing import List, Optional, Tuple

from flask import Response, jsonify, request
from flask_login import login_required
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook, sanitise_user_input

from . import novel_compound_bp


@novel_compound_bp.route("/_novel_compound", methods=["GET", "POST"])
@login_required
def novel_compound() -> Response:
    """
    Adds a novel compound to the database and links it to the workbook of the current reaction.
    The novel compound data is saved from the data the user enters in the form.
    This function is triggered either from the 'add new reagent/solvent to database'
    or when a novel structure is drawn in the sketcher.

    Returns:
        Flask Response with feedback about the operation.
    """
    # get the active workbook and verify user belongs
    workgroup_name = str(request.form["workgroup"])
    workbook_name = str(request.form["workbook"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)

    # validate the name is allowed and unique
    name, feedback = validate_name()
    if feedback:
        return jsonify({"feedback": feedback})
    if not is_name_unique(name, workbook):
        feedback = "A compound with this name is already in the database"
        return jsonify({"feedback": feedback})

    # validate numerical inputs
    expected_num_ls = validate_numbers()
    if expected_num_ls is None:
        return jsonify(
            {
                "feedback": "Molecular weight, density, and concentration must be empty or a positive number"
            }
        )
    density, concentration, mol_weight = expected_num_ls

    # validate the cas and check is unique
    cas_feedback = validate_cas(workbook)
    if cas_feedback:
        return jsonify({"feedback": cas_feedback})

    # If we have a SMILES, we use it to get additional representations.
    mol_formula, inchi, inchi_key = calculate_molecule_identifiers()

    # validate the hazards
    hazards = validate_hazards()

    # create novel compound
    nc = services.novel_compound.add(
        name,
        cas_feedback,
        mol_formula,
        mol_weight,
        density,
        concentration,
        hazards,
        request.form["smiles"],
        inchi,
        inchi_key,
        workbook.id,
    )

    # if the user has added the compound as a novel solvent, we add this to the solvent table too
    component_type = request.form["component"]
    if component_type == "solvent":
        services.solvent.add(name, hazards, nc)

    feedback = "Compound added to the database"
    return jsonify({"feedback": feedback})


def validate_name() -> Tuple[str, str]:
    """
    Validates and sanitizes the compound name from the request.

    Returns:
        Validated name or error message if validation fails.
    """
    name = sanitise_user_input(request.form["name"])
    feedback = None
    if len(name) > 200:
        feedback = "Name must be under 200 characters long"
    if not name:
        feedback = "Compound requires a name"
    return name, feedback


def is_name_unique(name: str, workbook: models.WorkBook) -> bool:
    """
    Checks if the compound name is unique within the given workbook.

    Args:
        name: The compound name to check.
        workbook: Workbook model.

    Returns:
        True if the name is unique, False otherwise.
    """
    name_check = services.novel_compound.get_novel_compound_from_name_and_workbook(
        name, workbook
    )
    second_name_check = services.compound.get_compound_from_name(name)
    return name_check is None and second_name_check is None


def validate_numbers() -> Optional[List[float]]:
    """
    Validates and extracts numerical values from the request.

    Returns:
        List of validated numerical values or None if validation fails.
    """
    expected_num_ls = [
        request.form["density"],
        request.form["concentration"],
        request.form["molWeight"],
    ]
    expected_num_ls = [float(x) if x != "" else None for x in expected_num_ls]
    for entry in expected_num_ls:
        if entry is not None and not check_positive_number(entry):
            return None
    return expected_num_ls


def validate_cas(workbook: models.WorkBook) -> Optional[str]:
    """
    Validates the CAS number from the request. Checks it fits the expected CAS pattern with a regular expression
    and checks it is unique within the compound database and the workbook's novel compound collection

    Args:
        workbook: Workbook we are checking to see if a compound with this cas is present in

    Returns:
        Error message if validation fails, None otherwise.
    """
    cas = sanitise_user_input(request.form["cas"])
    if cas:
        cas_regex = r"^[0-9]{1,7}-\d{2}-\d$"
        if not re.findall(cas_regex, cas):
            return "CAS invalid."
        cas_duplicate_check1 = services.compound.get_compound_from_cas(cas)
        cas_duplicate_check2 = (
            services.novel_compound.get_novel_compound_from_cas_and_workbook(
                cas, workbook
            )
        )
        if cas_duplicate_check1 or cas_duplicate_check2:
            return (
                "CAS already in database. Please add this compound to the reaction table"
                " by searching for the CAS in the reagent box"
            )
    return None


def check_positive_number(s: float) -> bool:
    """Checks the entry is a positive number."""
    try:
        return s >= 0
    except ValueError:
        return False


def calculate_molecule_identifiers() -> Tuple[str, Optional[str], Optional[str]]:
    """
    Calculate additional molecule identifiers if SMILES is present.

    Returns:
        Tuple containing molecule formula, InChI notation, and InChIKey.
    """
    smiles = request.form["smiles"]
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "", None, None
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        inchi = Chem.MolToInchi(mol)
        inchi_key = Chem.MolToInchiKey(mol)
    else:
        mol_formula = ""
        inchi = None
        inchi_key = None
    return mol_formula, inchi, inchi_key


def validate_hazards() -> str:
    """
    Check hazard codes are in the correct format.

    Returns:
        Hazard codes or "Unknown" if not provided or feedback for error message if invalid
    """
    hazards = sanitise_user_input(request.form["hPhrase"])
    if hazards:
        hazard_codes_ls = hazards.split("-")
        for hazard_code in hazard_codes_ls:
            hazard_match = services.hazard_code.get(hazard_code)
            if hazard_match is None:
                feedback = f'Hazard code "{hazard_code}" is invalid. Must be valid hazard code and formatted correctly. e.g., H200-H301.'
                return feedback
    else:
        hazards = "Unknown"
    return hazards
