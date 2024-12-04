from flask import (
    Response,
    abort,
    escape,
    flash,
    jsonify,
    redirect,
    render_template,
    request,
    url_for,
)
from flask_api import status
from flask_login import current_user, login_required
from sources import models, services
from sources.auxiliary import (
    abort_if_user_not_in_workbook,
    security_member_workgroup_workbook,
)
from sources.extensions import db
from sources.decorators import workbook_member_required

from . import reaction_list_bp


# delete reaction
@reaction_list_bp.route(
    "/delete_reaction/<reaction_id>/<workgroup>/<workbook>", methods=["GET", "POST"]
)
@login_required
@workbook_member_required
def delete_reaction(reaction_id: str, workgroup: str, workbook: str) -> Response:
    # must be logged in, a member of the workgroup and workbook and the creator of the reaction
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


@reaction_list_bp.route("/get_reactions", methods=["GET", "POST"])
@login_required
def get_reactions() -> Response:
    """Gets a list of reactions for the active workbook. Reaction data is sent as a list of dictionaries."""
    sort_crit = str(request.form["sortCriteria"])
    workbook_name = str(request.form["workbook"])
    workgroup_name = str(request.form["workgroup"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    reactions = services.reaction.list_active_in_workbook(
        workbook_name, workgroup_name, sort_crit
    )
    new_reaction_id = services.reaction.get_next_reaction_id_for_workbook(workbook.id)
    reactions = services.reaction.to_dict(reactions)
    reaction_details = render_template(
        "_saved_reactions.html",
        reactions=reactions,
        sort_crit=escape(sort_crit),
        workgroup=workgroup_name,
        new_reaction_id=new_reaction_id,
    )
    return jsonify({"reactionDetails": reaction_details})


@reaction_list_bp.route("/get_schemata", methods=["GET", "POST"])
@login_required
def get_schemata() -> Response:
    """
    Gets a list of reaction schemes for the active workbook. Images made from reaction SMILES using rdMolDraw2D (RDKit)
    """
    workbook_name = str(request.form["workbook"])
    workgroup_name = str(request.form["workgroup"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    size = str(request.form["size"])
    sort_crit = str(request.form["sortCriteria"])
    reaction_list = services.reaction.list_active_in_workbook(
        workbook_name, workgroup_name, sort_crit
    )
    schemes = services.reaction.make_scheme_list(reaction_list, size)
    return {"schemes": schemes, "sort_crit": escape(sort_crit)}


@reaction_list_bp.route("/get_new_reaction_id", methods=["GET", "POST"])
@login_required
def get_new_reaction_id() -> Response:
    workbook_name = request.json.get("workbook")
    workgroup_name = request.json.get("workgroup")

    if workgroup_name is None or workbook_name is None:
        return jsonify("Bad Request"), status.HTTP_400_BAD_REQUEST

    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )

    new_id = services.reaction.get_next_reaction_id_for_workbook(workbook.id)

    return jsonify(new_id)
