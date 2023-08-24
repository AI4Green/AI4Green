from flask import (Response, jsonify,
                   render_template, request, flash, redirect, url_for)
from flask_login import \
    login_required, current_user
from sources.extensions import db
from sources import models
from sources.auxiliary import security_member_workgroup_workbook

from . import reaction_list_bp  # imports the blueprint of route
from .reaction_list import get_reaction_list, get_scheme_list


# delete reaction
@reaction_list_bp.route(
    "/delete_reaction/<reaction_id>/<workgroup>/<workbook>", methods=["GET", "POST"]
)
@login_required
def delete_reaction(reaction_id: str, workgroup: str, workbook: str) -> Response:
    # must be logged in a member of the workgroup and workbook
    if not security_member_workgroup_workbook(workgroup, workbook):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    # find reaction
    reaction = (
        db.session.query(models.Reaction)
        .join(models.WorkBook)
        .join(models.WorkGroup)
        .join(models.Person)
        .join(models.User)
        .filter(models.Reaction.reaction_id == reaction_id)
        .filter(models.WorkBook.name == workbook)
        .filter(models.WorkGroup.name == workgroup)
        .filter(models.User.email == current_user.email)
        .first()
    )
    # check user is creator of reaction
    if not reaction:
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    # change to inactive
    reaction.status = "inactive"
    db.session.commit()
    return redirect(
        url_for(
            "workgroup.workgroup",
            workgroup_selected=workgroup,
            workbook_selected=workbook,
        )
    )


# Get reactions
@reaction_list_bp.route("/get_reactions", methods=["GET", "POST"])
@login_required
def get_reactions() -> Response:
    # must be logged in
    sort_crit = str(request.form["sortCriteria"])
    workbook = str(request.form["workbook"])
    workgroup = str(request.form["workgroup"])
    reactions = get_reaction_list(workbook, workgroup, sort_crit)
    reaction_details = render_template(
        "_saved_reactions.html", reactions=reactions, sort_crit=sort_crit
    )
    return jsonify({"reactionDetails": reaction_details})


@reaction_list_bp.route("/get_schemata", methods=["GET", "POST"])
@login_required
def get_schemata() -> Response:
    # must be logged in
    workbook = str(request.form["workbook"])
    workgroup = str(request.form["workgroup"])
    size = str(request.form["size"])
    sort_crit = str(request.form["sortCriteria"])
    schemes = get_scheme_list(workbook, workgroup, sort_crit, size)
    return {"schemes": schemes, "sort_crit": sort_crit}
