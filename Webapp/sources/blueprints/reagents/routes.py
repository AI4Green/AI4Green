import re

from flask import Response, jsonify, request
from sources import db, models, services
from sources.auxiliary import abort_if_user_not_in_workbook, sanitise_user_input
from sqlalchemy import func

# request parses incoming request data and gives access to it
# jsonify is used to send a JSON response to the browser
from . import reagents_bp  # imports the blueprint of the reagents route


# Getting the reagent name from browser and returning its data
@reagents_bp.route("/_reagents", methods=["POST"])
def reagents() -> Response:
    # must be logged in
    """This function gets a reagent name from browser,
    makes request to the reagent database, and returns
    its data back to show it in the reaction table"""
    reagent = sanitise_user_input(
        request.form["reagent"]
    )  # gets the reagent from browser
    reagent = reagent.replace("\u00A0", " ")  # add back in space to search database
    workgroup_name = str(request.form["workgroup"])
    workbook_name = str(request.form["workbook"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)

    novel_compound = False
    number = request.form["number"]  # gets the reagent number from the browser
    cas_regex = r"\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b"
    cas_number = re.findall(cas_regex, reagent)
    if cas_number:
        reagent = cas_number[0]
        found_reagent = (
            db.session.query(models.Compound)
            .filter(models.Compound.cas == reagent)
            .first()
        )
        if (
            found_reagent is None and workbook
        ):  # check the novel compound db for the cas number
            found_reagent = (
                db.session.query(models.NovelCompound)
                .filter(models.NovelCompound.cas == reagent)
                .join(models.WorkBook)
                .filter(models.WorkBook.id == workbook.id)
                .first()
            )
        if found_reagent is None:
            # if it's a cas but not in database then return for adding novel compound redirect
            return jsonify({"cas_not_found": True, "reagent": reagent})
    else:  # non-cas entry. We search reagent common names
        found_reagent = (
            db.session.query(models.Compound)
            .filter(func.lower(models.Compound.name) == reagent.lower())
            .first()
        )
        if found_reagent is None and workbook:  # check novel compound db
            found_reagent = (
                db.session.query(models.NovelCompound)
                .filter(func.lower(models.NovelCompound.name) == reagent.lower())
                .join(models.WorkBook)
                .filter(models.WorkBook.id == workbook.id)
                .first()
            )
            if found_reagent:
                novel_compound = True
        if (
            found_reagent is None
        ):  # If still no match find a list of partial matches and send to the frontend
            reagent_list = (
                db.session.query(models.Compound)
                .filter(func.lower(models.Compound.name).startswith(reagent.lower()))
                .order_by(models.Compound.name.asc())
                .limit(100)
                .all()
            )
            reagent_names = [x.name for x in reagent_list]
            return jsonify(
                {
                    "identifiers": reagent_names,
                    "reagent_not_found": True,
                    "reagent": reagent,
                }
            )
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
        if novel_compound is True:
            primary_key = f"('{found_reagent.name}', {found_reagent.workbook})"
        else:
            primary_key = found_reagent.id

        smiles = found_reagent.smiles
        # sends the JSON with the reagent data back to the browser
        return jsonify(
            {
                "name": reagent_name,
                "molWeight": mol_weight,
                "number": number,
                "density": density,
                "concentration": concentration,
                "hazards": hazards,
                "primary_key": primary_key,
                "smiles": smiles,
            }
        )
    return "", 204
