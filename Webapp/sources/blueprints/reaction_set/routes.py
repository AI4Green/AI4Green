from typing import Dict, List, Union

import requests
from flask import current_app, jsonify, redirect, render_template, request, url_for
from flask_login import current_user
from sources import models, services

from . import reaction_set_bp


@reaction_set_bp.route("/reaction_set/<workgroup_name>/<workbook_name>/<set_name>")
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
