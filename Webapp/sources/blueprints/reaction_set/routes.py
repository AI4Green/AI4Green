from typing import Dict, List, Union

from flask import current_app, jsonify, redirect, render_template, request, url_for
from flask_login import current_user, login_required
from sources import models, services

from . import reaction_set_bp


@reaction_set_bp.route("/reaction_set/<workgroup_name>/<workbook_name>/<set_name>")
@login_required
def reaction_set(set_name, workgroup_name, workbook_name):
    # to do
    # add workbook, workgroup to reaction set page
    # fix atuosave sketcher to only update reaction table in set mode (do we want to try autosaving?)
    # then fix apply to well and apply to all
    # colours for unsaved/edited wells
    #

    r_set = services.reaction_set.get_from_names(
        set_name, workgroup_name, workbook_name
    )

    serialised_set = serialise_reaction_set(r_set)

    return render_template(
        "reaction_set.html",
        reaction_set=serialised_set,
        number_of_reactions=len(r_set.reactions),
        workgroup_name=workgroup_name,
        workbook_name=workbook_name,
        reactor_dimensions=r_set.reactor_dimensions,
    )


def serialise_reaction_set(reaction_set: models.ReactionSet):
    return {
        "id": reaction_set.id,
        "name": reaction_set.name,
        "reactions": [serialise_reaction(r) for r in reaction_set.reactions],
    }


def serialise_reaction(reaction: models.Reaction) -> Dict[str, str]:
    return {
        "reaction_id": reaction.reaction_id,
        "name": reaction.name,
        "smiles": reaction.reaction_smiles,
        # maybe more
    }


@reaction_set_bp.route("/click_and_drag")
def click_and_drag():
    return render_template("click-and-drag.html")


@reaction_set_bp.route("/new_reaction_set", methods=["GET", "POST"])
@login_required
def new_reaction_set():
    data = request.get_json()
    workgroup_name = data.get("workgroup", None)
    workbook_name = data.get("workbook", None)
    reactor_dimensions = data.get("reactorDimensions", None)
    # probs need some checks here to make sure necessary data is present

    new_set_id = data.get("setID", None)
    new_set_name = data.get("setName")

    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    creator = services.person.from_current_user_email()

    number_of_reactions = get_number_of_reactions(reactor_dimensions)

    # new_set_id = services.reaction_set.next_id_in_workbook(workbook.id)
    reaction_table = services.reaction.empty_reaction_table()
    summary_table = services.summary.empty_summary_table()

    reactions = []
    for i in range(number_of_reactions):
        # reaction_id = services.reaction.get_next_reaction_id_for_workbook(workbook.id)
        reactions.append(
            services.reaction.add(
                name=new_set_name + "-" + str(i + 1),
                reaction_id=new_set_id + "-" + str(i + 1),
                creator=current_user,
                workbook_id=workbook.id,
                reaction_table=reaction_table,
                summary_table=summary_table,
                reaction_smiles="",
            )
        )

    services.reaction_set.add(
        name=new_set_name,
        set_id=new_set_id,
        creator=creator,
        workbook=workbook,
        reactions=reactions,
        reactor_dimensions=reactor_dimensions,
    )

    return jsonify({"feedback": "Set Created!"})


def get_number_of_reactions(reactor_dimensions):
    if reactor_dimensions.get("reactorType") == "carousel":
        return reactor_dimensions.get("numberOfReactions")
    elif reactor_dimensions.get("reactorType") == "well-plate":
        return reactor_dimensions.get("rows") * reactor_dimensions.get("columns")
