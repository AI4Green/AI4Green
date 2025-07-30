"""
This module receives a reaction from Marvin JS as a
GET request and renders the reaction table template
"""

from flask import json, jsonify, render_template, request
from flask_login import login_required
from sources import models, services
from sources.dto import ReactionNoteSchema

from . import reaction_table_bp


@reaction_table_bp.route("/autoupdate_reaction_table", methods=["GET", "POST"])
# @workbook_member_required
def autoupdate_reaction_table():
    """
    I guess the idea is to take any smiles from the front end and generate the reaction table

    Needs to get reaction from db too, but ignore for now
    """
    # get user workbook
    demo = request.json.get("demo")
    tutorial = request.json.get("tutorial")
    workbook = None
    workgroup = None
    reaction = None
    polymer_indices = []
    if demo != "demo" and tutorial != "yes":
        workgroup_name = request.json.get("workgroup")
        workbook_name = request.json.get("workbook")
        reaction_id = request.json.get("reaction_id")

        workgroup = services.workgroup.from_name(workgroup_name)
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            workgroup_name, workbook_name
        )
        reaction = services.reaction.get_from_reaction_id_and_workbook_id(
            reaction_id, workbook.id
        )

    reaction_smiles = request.json.get("reaction_smiles")
    if not reaction_smiles:
        return jsonify({"error": "Missing data!"})

    polymer_indices = request.json.get("polymer_indices")

    reaction_table = services.reaction_table.ReactionTable(
        reaction, workgroup, workbook, demo, tutorial
    )

    reaction_table.update(reaction_smiles, polymer_indices)

    return reaction_table.render()


@reaction_table_bp.route("/reload_reaction_table", methods=["GET", "POST"])
def reload_reaction_table():
    # load variables from front end
    workbook = request.json.get("workbook")
    workgroup = request.json.get("workgroup")
    reaction_id = request.json.get("reaction_id")
    demo = request.json.get("demo")
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup, workbook
    )
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook.id
    )

    reaction_table = services.reaction_table.ReactionTable(
        reaction, workgroup, workbook, demo, "no"
    )
    print(reaction_table.reactants, reaction_table.products)

    return reaction_table.render()


@reaction_table_bp.route("/_save_reaction_note", methods=["POST"])
@login_required
@reaction_table_bp.doc(security="sessionAuth")
def save_reaction_note():
    """
    Saves an reaction_note to the reaction object

    Returns:
        flask.Response: A JSON response with the reaction_note object
    """
    workgroup_name = request.form["workgroup"]
    workbook_name = request.form["workbook"]
    workbook_object = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    reaction_id = request.form["reactionID"]
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook_object.id
    )
    reaction_note_text = request.form["reactionNoteText"]
    author = services.person.from_current_user_email()
    new_addendum = services.reaction.add_addendum(reaction, reaction_note_text, author)
    schema = ReactionNoteSchema()
    return jsonify({"reaction_note": schema.dump(new_addendum)})
