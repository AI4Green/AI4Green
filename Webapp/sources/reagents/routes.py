from flask import request, jsonify
# request parses incoming request data and gives access to it
# jsonify is used to send a JSON response to the browser
from sources.reagents import reagents_bp  # imports the blueprint of the reagents route
from sources import auxiliary, db  # imports the module with auxiliary functions
from pony.orm import select
import re
from flask_login import login_required


# Getting the reagent name from browser and returning its data
@reagents_bp.route('/_reagents', methods=['POST'])
@login_required
def reagents():
    # must be logged in
    """This function gets a reagent name from browser,
    makes request to the reagent database, and returns
    its data back to show it in the reaction table"""
    reagent = auxiliary.sanitise_user_input(request.form['reagent'])  # gets the reagent from browser
    reagent = reagent.replace("\u00A0", " ")  # add back in space to search database
    workgroup = str(request.form['workgroup'])
    workbook_name = str(request.form['workbook'])
    workbook = select(b for b in db.WorkBook if b.name == workbook_name and b.group.name == workgroup).first()
    number = request.form['number']  # gets the reagent number from the browser
    cas_regex = r'\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b'
    if re.findall(cas_regex, reagent):  # if it's a cas, find the corresponding name in the compound db
        found_reagent = select(c for c in db.Compound if reagent == c.cas).first()
        if found_reagent is None: # check the novel compound db for the cas number
            found_reagent = select(c for c in db.NovelCompound if reagent == c.cas and c.workbook == workbook).first()
        if found_reagent is None:
            # if it's a cas but not in database then return for adding novel compound redirect
            return jsonify({'cas_not_found': True, 'reagent': reagent})
    else:  # non-cas entry. We search reagent common names
        found_reagent = select(c for c in db.Compound if reagent.lower() == c.name.lower()).first()
        if found_reagent is None:  # check novel compound db
            found_reagent = select(c for c in db.NovelCompound if reagent.lower() == c.name.lower() and c.workbook == workbook).first()
        if found_reagent is None:  # If still no match find a list of partial matches and send to the frontend
            a = select(c.name for c in db.Compound if c.name.lower().startswith(reagent.lower()))[:100]
            c = sorted(list(set(a)))[:100]
            return jsonify({'identifiers': c, 'reagent_not_found': True, 'reagent': reagent})
    if found_reagent is not None:
        reagent_name = found_reagent.name  # assign reagent to name
        # then its hazard phrase, density, molar weight,
        # and concentration are selected from the reagent
        mol_weight = found_reagent.molec_weight  # reagent molar weight
        density = found_reagent.density  # reagent density
        density = 0 if density == "-" else density
        # assigns 0 if the reagent density is not found in the database
        concentration = found_reagent.concentration  # reagent concentration
        concentration = 0 if concentration == "-" else concentration
        # assigns 0 if the reagent density is not found in the database
        hazards = found_reagent.hphrase  # reagent hazards
        primary_key = found_reagent.id
        smiles = found_reagent.smiles
        # sends the JSON with the reagent data back to the browser
        return jsonify({'name': reagent_name, 'molWeight': mol_weight, 'number': number, 'density': density,
                        'concentration': concentration, 'hazards': hazards, 'primary_key': primary_key,
                        'smiles': smiles})
    return '', 204
