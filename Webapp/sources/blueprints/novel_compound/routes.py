import re
from datetime import datetime

import pytz
from flask import Response, jsonify, request
from flask_login import current_user, login_required
from rdkit import Chem  # Used for converting smiles to inchi
from rdkit.Chem import rdMolDescriptors
from sqlalchemy import func

# render_template renders html templates
# request parses incoming request data and gives access to it
# jsonify is used to send a JSON response to the browser
from sources import (auxiliary,  # imports the module with auxiliary functions
                     db, models)

from . import \
    novel_compound_bp  # imports the blueprint of the reaction table route


@novel_compound_bp.route("/_novel_compound", methods=["GET", "POST"])
@login_required
def novel_compound() -> Response:
    # must be logged in
    # get user, their workbook, and the chemicals within that workbook to check the name is unique
    name = auxiliary.sanitise_user_input(request.form["name"])
    if not name:
        feedback = "Compound requires a name"
        return jsonify({"feedback": feedback})
    workgroup = str(request.form["workgroup"])
    workbook_name = str(request.form["workbook"])
    workbook = (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.name == workbook_name)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )
    # check novel compound db
    name_check = (
        db.session.query(models.NovelCompound)
        .filter(func.lower(models.NovelCompound.name) == name.lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .first()
    )
    # check compound db
    second_name_check = (
        db.session.query(models.Compound)
        .filter(func.lower(models.Compound.name) == name.lower())
        .first()
    )
    # name must be unique within workbook
    if name_check or second_name_check:
        feedback = "A compound with this name is already in the database"
        return jsonify({"feedback": feedback})

    # if values are provided mol_weight, density, and conc must be >= 0
    density = request.form["density"]
    concentration = request.form["concentration"]
    mol_weight = request.form["molWeight"]
    expected_num_ls = [density, concentration, mol_weight]
    # turning empty strings into None to fit database constraints
    expected_num_ls = [x if x != "" else None for x in expected_num_ls]
    for entry in expected_num_ls:
        # if not empty or 0, it must be checked.
        if entry is not None:
            valid = check_positive_number(entry)
            if valid is False:
                feedback = "Molecular weight, density, and concentration must be empty or a positive number"
                return jsonify({"feedback": feedback})
    # unpack list
    density, concentration, mol_weight = expected_num_ls
    # if cas provided, must be valid
    cas = auxiliary.sanitise_user_input(request.form["cas"])
    if cas:
        cas_regex = r"^[0-9]{1,7}-\d{2}-\d$"
        if not re.findall(cas_regex, cas):
            feedback = "CAS invalid."
            return jsonify({"feedback": feedback})
        cas_duplicate_check1 = (
            db.session.query(models.Compound).filter(models.Compound.cas == cas).first()
        )
        cas_duplicate_check2 = (
            db.session.query(models.NovelCompound)
            .filter(models.Compound.cas == cas)
            .join(models.WorkBook)
            .filter(models.WorkBook.id == workbook.id)
            .first()
        )
        if cas_duplicate_check1 or cas_duplicate_check2:
            return jsonify(
                {
                    "feedback": "CAS already in database. Please add this compound to the reaction table"
                    " by searching for the CAS in the reagent box"
                }
            )
    # calculate additional molecule identifiers if smiles is present
    smiles = request.form["smiles"]
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return jsonify({"feedback": "Invalid smiles"})
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        inchi = Chem.MolToInchi(mol)
        inchi_key = Chem.MolToInchiKey(mol)
    else:
        mol_formula = ""
        inchi = None
        inchi_key = None
    # check hazard codes are in correct format
    hazards = auxiliary.sanitise_user_input(request.form["hPhrase"])
    if hazards:
        hazards_ls = hazards.split("-")
        for hazard in hazards_ls:
            hazard_match = (
                db.session.query(models.HazardCode)
                .filter(models.HazardCode.code == hazard)
                .first()
            )
            if hazard_match is None:
                feedback = f'Hazard code "{hazard}" is invalid. Must be valid hazard code and formatted correctly. e.g., H200-H301.'
                return jsonify({"feedback": feedback})
    else:
        hazards = "Unknown"
    current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)

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
        workbook=workbook.id,
        time_of_creation=current_time,
    )
    db.session.add(nc)
    db.session.commit()
    component_type = request.form["component"]
    if component_type == "solvent":
        model = models.Solvent(
            name=name,
            flag=5,
            hazard=hazards,
            novel_compound=[nc],
            time_of_creation=current_time,
        )
        db.session.add(model)
        db.session.commit()
    feedback = "Compound added to the database"
    return jsonify({"feedback": feedback})


def check_positive_number(s: str) -> bool:
    """Checks the entry is a positive number"""
    try:
        return float(s) >= 0
    except ValueError:
        return False
