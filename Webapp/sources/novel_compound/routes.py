from flask import request, jsonify
from flask_login import login_required, current_user
# render_template renders html templates
# request parses incoming request data and gives access to it
# jsonify is used to send a JSON response to the browser
from sources import auxiliary, db  # imports the module with auxiliary functions
from sources.novel_compound import novel_compound_bp  # imports the blueprint of the reaction table route
from rdkit import Chem  # Used for converting smiles to inchi
from rdkit.Chem import rdMolDescriptors
from pony.orm import select
import re


@novel_compound_bp.route('/_novel_compound', methods=['GET', 'POST'])
@login_required
def novel_compound():
    # must be logged in
    # get user, their workbook, and the chemicals within that workbook to check the name is unique
    name = auxiliary.sanitise_user_input(request.form['name'])
    if not name:
        feedback = 'Compound requires a name'
        return jsonify({'feedback': feedback})
    workgroup = str(request.form['workgroup'])
    workbook_name = str(request.form['workbook'])
    workbook = select(b for b in db.WorkBook if b.name == workbook_name and b.group.name == workgroup).first()
    # check novel compound db
    name_check = select(x for x in db.NovelCompound if x.workbook == workbook and x.name.lower() == name.lower()).first()
    # check compound db
    second_name_check = select(x for x in db.Compound if x.name.lower() == name.lower()).first()
    # name must be unique within workbook
    if name_check or second_name_check:
        feedback = 'A compound with this name is already in the database'
        return jsonify({'feedback': feedback})

    # if values are provided mol_weight, density, and conc must be >= 0
    density = request.form['density']
    concentration = request.form['concentration']
    mol_weight = request.form['molWeight']
    expected_num_ls = [density, concentration, mol_weight]
    # turning empty strings into None to fit database constraints
    expected_num_ls = [x if x != '' else None for x in expected_num_ls]
    for entry in expected_num_ls:
        # if not empty or 0, it must be checked.
        if entry is not None:
            valid = check_positive_number(entry)
            if valid is False:
                feedback = 'Molecular weight, density, and concentration must be empty or a positive number'
                return jsonify({'feedback': feedback})
    # unpack list
    density, concentration, mol_weight = expected_num_ls
    # if cas provided, must be valid
    cas = auxiliary.sanitise_user_input(request.form['cas'])
    if cas:
        cas_regex = r'^[0-9]{1,7}-\d{2}-\d$'
        if not re.findall(cas_regex, cas):
            feedback = 'CAS invalid.'
            return jsonify({'feedback': feedback})
        cas_duplicate_check1 = select(x for x in db.Compound if x.cas == cas).first()
        cas_duplicate_check2 = select(x for x in db.NovelCompound if x.cas == cas and x.workbook == workbook).first()
        if cas_duplicate_check1 or cas_duplicate_check2:
            return jsonify({'feedback': 'CAS already in database. Please add this compound to the reaction table'
                                        ' by searching for the CAS in the reagent box'})
    # calculate additional molecule identifiers if smiles is present
    smiles = request.form['smiles']
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return jsonify({'feedback': 'Invalid smiles'})
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        inchi = Chem.MolToInchi(mol)
        inchi_key = Chem.MolToInchiKey(mol)
    else:
        mol_formula = ''
        inchi = None
        inchi_key = None
    # check hazard codes are in correct format
    hazards = auxiliary.sanitise_user_input(request.form['hPhrase'])
    if hazards:
        hazards_ls = hazards.split('-')
        for hazard in hazards_ls:
            hazard_match = select(x for x in db.HazardCode if x.code == hazard).first()
            if hazard_match is None:
                feedback = f'Hazard code "{hazard}" is invalid. Must be valid hazard code and formatted correctly. e.g., H200-H301.'
                return jsonify({'feedback': feedback})
    else:
        hazards = 'Unknown'

    nc = db.NovelCompound(name=name, cas=cas, molec_formula=mol_formula, molec_weight=mol_weight, density=density,
                          concentration=concentration,
                          hphrase=hazards, smiles=smiles, InChI=inchi, InChIKey=inchi_key, workbook=workbook)
    component_type = request.form['component']
    if component_type == 'solvent':
        db.Solvent(name=name, flag=5, hazard=hazards, compound=None, novel_compound=nc)
    feedback = 'Compound added to the database'
    return jsonify({'feedback': feedback})


def check_positive_number(s):
    """Checks the entry is a positive number"""
    try:
        if float(s) >= 0:
            return True
        else:
            return False
    except ValueError:
        return False
