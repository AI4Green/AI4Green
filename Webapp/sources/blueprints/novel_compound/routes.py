import re
from datetime import datetime
from typing import List, Optional, Tuple, Union

import pytz
from flask import Response, jsonify, request
from flask_login import current_user, login_required
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from sources import db, models
from sources.auxiliary import (
    abort_if_user_not_in_workbook,
    get_workbook_from_group_book_name_combination,
    sanitise_user_input,
)
from sqlalchemy import func

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
    name, feedback = validate_name(request)
    if feedback:
        return jsonify({"feedback": feedback})

    workgroup_name = str(request.form["workgroup"])
    workbook_name = str(request.form["workbook"])
    workbook = get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)

    if not is_name_unique(name, workbook):
        feedback = "A compound with this name is already in the database"
        return jsonify({"feedback": feedback})

    expected_num_ls = validate_numbers(request)
    if expected_num_ls is None:
        return jsonify(
            {
                "feedback": "Molecular weight, density, and concentration must be empty or a positive number"
            }
        )

    density, concentration, mol_weight = expected_num_ls

    cas_feedback = validate_cas(request, workbook)
    if cas_feedback:
        return jsonify({"feedback": cas_feedback})

    mol_formula, inchi, inchi_key = calculate_molecule_identifiers(request)

    hazards = validate_hazards(request)

    current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)

    nc = create_novel_compound(
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
        current_time,
    )

    component_type = request.form["component"]
    if component_type == "solvent":
        create_solvent_model(name, hazards, nc, current_time)

    feedback = "Compound added to the database"
    return jsonify({"feedback": feedback})


def validate_name(request: request) -> Tuple[Optional[str], Optional[str]]:
    """
    Validates and sanitizes the compound name from the request.

    Args:
        request: Flask request object.

    Returns:
        Tuple of (validated_name, feedback) where feedback is None if validation is successful.
    """
    name = sanitise_user_input(request.form["name"])
    if len(name) > 200:
        return None, "Name must be under 200 characters long"
    if not name:
        return None, "Compound requires a name"
    return name, None


def is_name_unique(name: str, workbook: models.WorkBook) -> bool:
    """
    Checks if the compound name is unique within the specified workbook.

    Args:
        name: Compound name.
        workbook: Workbook model.

    Returns:
        True if the name is unique, False otherwise.
    """
    name_check = get_novel_compound_by_name_and_workbook(name, workbook)
    second_name_check = get_compound_by_name(name)
    return not (name_check or second_name_check)


def validate_numbers(request: request) -> Union[List[Optional[float]]]:
    """
    Validates and extracts numerical values from the request.

    Args:
        request: Flask request object.

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
            return
    return expected_num_ls


def validate_cas(request: request, workbook: models.WorkBook) -> Optional[str]:
    """
    Validates the CAS number from the request.

    Args:
        request: Flask request object.
        workbook: Workbook model.

    Returns:
        Error message if validation fails, None otherwise.
    """
    cas = sanitise_user_input(request.form["cas"])
    if cas:
        cas_regex = r"^[0-9]{1,7}-\d{2}-\d$"
        if not re.findall(cas_regex, cas):
            return "CAS invalid."
        cas_duplicate_check1 = get_compound_by_cas(cas)
        cas_duplicate_check2 = get_novel_compound_by_cas_and_workbook(cas, workbook)
        if cas_duplicate_check1 or cas_duplicate_check2:
            return (
                "CAS already in database. Please add this compound to the reaction table"
                " by searching for the CAS in the reagent box"
            )
    return None


def calculate_molecule_identifiers(
    request: request,
) -> Tuple[str, Optional[str], Optional[str]]:
    """
    Calculates molecule identifiers from the SMILES provided in the request.

    Args:
        request: Flask request object.

    Returns:
        Tuple of molecule formula, InChI, and InChIKey.
    """
    smiles = request.form["smiles"]
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "Invalid smiles", None, None
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        inchi = Chem.MolToInchi(mol)
        inchi_key = Chem.MolToInchiKey(mol)
    else:
        mol_formula, inchi, inchi_key = "", None, None
    return mol_formula, inchi, inchi_key


def validate_hazards(request: request) -> str:
    """
    Validates hazard codes from the request.

    Args:
        request: Flask request object.

    Returns:
        Validated hazard codes or a default value if hazards are not provided.
    """
    hazards = sanitise_user_input(request.form["hPhrase"])
    if hazards:
        hazards_ls = hazards.split("-")
        for hazard in hazards_ls:
            hazard_match = (
                db.session.query(models.HazardCode)
                .filter(models.HazardCode.code == hazard)
                .first()
            )
            if hazard_match is None:
                return f'Hazard code "{hazard}" is invalid. Must be a valid hazard code and formatted correctly. e.g., H200-H301.'
    return "Unknown"


def create_novel_compound(
    name: str,
    cas: Optional[str],
    mol_formula: str,
    mol_weight: float,
    density: float,
    concentration: float,
    hazards: str,
    smiles: str,
    inchi: Optional[str],
    inchi_key: Optional[str],
    workbook_id: int,
    current_time: datetime,
) -> models.NovelCompound:
    """
    Creates a novel compound in the database.

    Args:
        name: Compound name.
        cas: CAS number.
        mol_formula: Molecule formula.
        mol_weight: Molecular weight.
        density: Density.
        concentration: Concentration.
        hazards: Hazard codes.
        smiles: SMILES notation.
        inchi: InChI notation.
        inchi_key: InChIKey.
        workbook_id: Workbook ID.
        current_time: Current timestamp.

    Returns:
        NovelCompound model.
    """
    nc = models.NovelCompound(
        name=name,
        cas=cas,
        molec_formula=mol_formula,
        molec_weight=mol_weight,
        density=density,
        concentration=concentration,
        hphrase=hazards,
        smiles=smiles,
        inchi=inchi,
        inchikey=inchi_key,
        workbook=workbook_id,
        time_of_creation=current_time,
    )
    db.session.add(nc)
    db.session.commit()
    return nc


def create_solvent_model(
    name: str, hazards: str, nc: models.NovelCompound, current_time: datetime
) -> models.Solvent:
    """
    Creates a solvent model in the database.

    Args:
        name: Solvent name.
        hazards: Hazard codes.
        nc: NovelCompound model.
        current_time: Current timestamp.

    Returns:
        Solvent model.
    """
    model = models.Solvent(
        name=name,
        flag=5,
        hazard=hazards,
        novel_compound=[nc],
        time_of_creation=current_time,
    )
    db.session.add(model)
    db.session.commit()
    return model


def check_positive_number(s: str) -> bool:
    """
    Checks if the entry is a positive number.

    Args:
        s: Input string.

    Returns:
        True if the entry is a positive number, False otherwise.
    """
    try:
        return float(s) >= 0
    except ValueError:
        return False


# SQLAlchemy queries
def get_novel_compound_by_name_and_workbook(
    name: str, workbook: models.WorkBook
) -> models.NovelCompound:
    return (
        db.session.query(models.NovelCompound)
        .filter(func.lower(models.NovelCompound.name) == name.lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )


def get_compound_by_name(name: str) -> models.Compound:
    return (
        db.session.query(models.Compound)
        .filter(func.lower(models.Compound.name) == name.lower())
        .first()
    )


def get_compound_by_cas(cas: str) -> models.Compound:
    return db.session.query(models.Compound).filter(models.Compound.cas == cas).first()


def get_novel_compound_by_cas_and_workbook(
    cas: str, workbook: models.WorkBook
) -> models.NovelCompound:
    return (
        db.session.query(models.NovelCompound)
        .filter(models.NovelCompound.cas == cas)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )
