from flask import request, jsonify
# request parses incoming request data and gives access to it
# jsonify is used to send a JSON response to the browser
from sources.solvents import solvents_bp  # imports the blueprint of the reagents route
from sources import db, auxiliary
from pony.orm import select
import re
from flask_login import login_required


# Getting the solvent name from browser and returning its data
@solvents_bp.route('/_solvents', methods=['POST'])
@login_required
def solvents():
    """This function gets a solvent name from browser,
    makes request to the solvent database, and returns
    its data back to show it in the reaction table"""
    solvent = auxiliary.sanitise_user_input(request.form['solvent'])  # gets the solvent from browser
    workgroup = request.form['workgroup']
    workbook_name = request.form['workbook']
    workbook = select(b for b in db.WorkBook if b.name == workbook_name and b.group.name == workgroup).first()
    number = request.form['number']  # gets the solvent number from the browser
    flag_rate = {1: 'hazard-highly-hazardous', 4: 'hazard-acceptable', 3: 'hazard-warning', 2: 'hazard-hazardous',
                 5: 'non-chem21'}  # flag rate dictionary
    cas_regex = r'\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b'
    flag_color = 'white'
    hazards = ''
    primary_key = ''
    alert_message = ''
    new_solvent = False
    if re.findall(cas_regex, solvent):  # if it's a cas, find the corresponding name
        solvent_cas_match = select(c for c in db.Compound if solvent == c.cas).first()
        if solvent_cas_match is None: # if not in the database then check novel compound db
            # alert_message = 'The solvent is not found in the PubChem database'
            novel_compound_match = select(c for c in db.NovelCompound if solvent == c.cas and c.workbook == workbook).first()
            if novel_compound_match is not None:
                primary_key = (novel_compound_match.name, str(novel_compound_match.workbook))
                hazards = novel_compound_match.hphrase
                solvent = novel_compound_match.name
                flag = novel_compound_match.solvent.flag  # solvent flag rate
                flag_color = flag_rate[flag]  # solvent flag colour

            if novel_compound_match is None:
                new_solvent = True
                alert_message = 'The solvent is not found in the available databases. You can add it as a new compound.'
                return jsonify({'num': number, 'solvent': solvent, 'flag': 'white', 'new_solvent': new_solvent,
                                'hazards': '', 'primary_key': primary_key, 'alert_message': alert_message})
        if solvent_cas_match is not None:
            primary_key = solvent_cas_match.id
            hazards = solvent_cas_match.hphrase  # solvent hazards
            solvent = solvent_cas_match.name  # assign solvent to name
            solvent_dropdown_match = solvent_cas_match.solvent  # check if linked to solvent in chem21 list
            if solvent_dropdown_match is not None:  # if solvent is in solvent db then we take solvent flag from there
                flag = solvent_dropdown_match.flag  # solvent flag rate
                flag_color = flag_rate[flag]  # solvent flag colour

    else:  # if not a cas string check if it's a solvent name from the dropdown list, then go directly to solvent db
        # look for solvent without novel compound i.e., chem21 default or novel compound within user's workbook
        solvent_dropdown_match = select(s for s in db.Solvent if solvent.lower() == s.name.lower() and
                                        (s.novel_compound is None or s.novel_compound.workbook == workbook)).first()
        if solvent_dropdown_match is not None:
            flag = solvent_dropdown_match.flag  # solvent flag rate
            flag_color = flag_rate[flag]  # solvent flag colour
            hazards = solvent_dropdown_match.hazard  # solvent hazards
            # if not novel compound then use compound id otherwise use novel compound composite key as pk
            if solvent_dropdown_match.novel_compound is None:
                # if no matching entry in compound db use 0
                if solvent_dropdown_match.compound:
                    primary_key = solvent_dropdown_match.compound.id
                else:
                    primary_key = 0
            else:
                primary_key = (solvent_dropdown_match.novel_compound.name, workbook.name)
        if solvent_dropdown_match is None:  # if not in the solvent dropdown then check the compound database for the name
            compound_db_match = select(c for c in db.Compound if solvent.lower() == c.name.lower()).first()
            if compound_db_match is not None:
                flag_color = 'non-chem21'
                hazards = compound_db_match.hphrase
                primary_key = compound_db_match.id
            if compound_db_match is None:  # if no match in the compound database then return nothing
                return '', 204
    # set pk to 0 if solvent does not have a pk.
    primary_key = 0 if primary_key is None else primary_key
    # sends the JSON with the solvent data back to the browser
    return jsonify({'num': number, 'solvent': solvent, 'flag': flag_color, 'hazards': hazards,
                    'primary_key': primary_key, 'alert_message': alert_message, 'new_solvent': new_solvent})
