from flask import (
    Response,
    abort,
    escape,
    flash,
    json,
    jsonify,
    redirect,
    render_template,
    request,
    url_for,
)
from flask_api import status
from flask_login import current_user, login_required
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook
from sources.decorators import workbook_member_required
from sources.extensions import db

from . import reaction_list_bp


# delete reaction
@reaction_list_bp.route(
    "/delete_reaction/<reaction_id>/<workgroup>/<workbook>", methods=["GET", "POST"]
)
@login_required
@workbook_member_required
@reaction_list_bp.doc(
    "sessionAuth",
    description="Requires login and membership in the specified workbook.",
)
def delete_reaction(reaction_id: str, workgroup: str, workbook: str) -> Response:
    """
    must be logged in, a member of the workgroup and workbook and the creator of the reaction
    Deletes a reaction from the list of reactions in the workbook by changing its status to 'inactive'.

    Args:
        reaction_id: the id of the reaction to be deleted
        workgroup: the name of the workgroup that the workbook belongs to
        workbook: the name of the workbook that the reaction belongs to
    Returns:
        flask.Response: Redirects to the workgroup page


    """
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
@reaction_list_bp.doc("sessionAuth")
def get_reactions() -> Response:
    """
    Gets a list of reactions for the active workbook. Reaction data is sent as a list of dictionaries.

    Returns:
        flask.Response: A JSON response with the reaction data
    """
    sort_crit = str(request.form.get("sortCriteria"))
    workbook_name = str(request.form.get("workbook"))
    workgroup_name = str(request.form.get("workgroup"))
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


@reaction_list_bp.route("/get_reaction_images", methods=["GET", "POST"])
@login_required
@workbook_member_required
@reaction_list_bp.doc(
    "sessionAuth",
    description="Requires login and membership in the specified workbook.",
)
def get_reaction_images(workgroup, workbook) -> Response:
    """
    Gets a list of reaction images for the active workbook.

    Returns:
        flask.Response: A JSON response with the reaction images
    """
    sort_crit = str(request.form.get("sortCriteria"))
    reaction_list = services.reaction.list_active_in_workbook(
        workbook, workgroup, sort_crit
    )
    reaction_images = services.reaction.make_reaction_image_list(reaction_list)
    return {"reaction_images": reaction_images, "sort_crit": escape(sort_crit)}


@reaction_list_bp.route("/get_smiles", methods=["GET", "POST"])
@login_required
@workbook_member_required
@reaction_list_bp.doc(
    "sessionAuth",
    description="Requires login and membership in the specified workbook.",
)
def get_smiles(workgroup, workbook) -> Response:
    """
    Gets a list of reaction smiles for the active workbook.

    Returns:
        flask.Response: A JSON response with the reaction smiles
    """
    sort_crit = str(request.form.get("sortCriteria"))
    reaction_list = services.reaction.list_active_in_workbook(
        workbook, workgroup, sort_crit
    )
    smiles_list = [r.reaction_smiles for r in reaction_list]
    return {"smiles": smiles_list, "sort_crit": escape(sort_crit)}


@reaction_list_bp.route("/get_rxns", methods=["GET", "POST"])
@login_required
@workbook_member_required
@reaction_list_bp.doc(
    "sessionAuth",
    description="Requires login and membership in the specified workbook.",
)
def get_rxns(workgroup, workbook) -> Response:
    """
    Gets a list of reaction RXNs for the active workbook.

    Returns:
        flask.Response: A JSON response with the reaction RXNs
    """
    sort_crit = str(request.form.get("sortCriteria"))
    reaction_list = services.reaction.list_active_in_workbook(
        workbook, workgroup, sort_crit
    )
    rxn_list = [r.reaction_rxn for r in reaction_list]
    return {"rxns": rxn_list, "sort_crit": escape(sort_crit)}


@reaction_list_bp.route("/_save_new_images", methods=["POST"])
@login_required
@reaction_list_bp.doc("sessionAuth")
def save_new_images():
    """
    Updates reaction dict with images

    Returns:
        flask.Response: A JSON response with feedback
    """
    # must be logged in
    workbook_name = str(request.form.get("workbook"))
    workgroup_name = str(request.form.get("workgroup"))
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    sort_crit = str(request.form.get("sortCriteria"))
    reaction_list = services.reaction.list_active_in_workbook(
        workbook_name, workgroup_name, sort_crit
    )
    images = json.loads(request.form.get("images"))
    for idx, reaction in enumerate(reaction_list):
        update_dict = {"reaction_image": images[idx]}
        reaction.update(**update_dict)
    feedback = "Reaction Updated!"
    return jsonify({"feedback": feedback})


@reaction_list_bp.route("/get_new_reaction_id", methods=["GET", "POST"])
@login_required
@reaction_list_bp.doc("sessionAuth")
def get_new_reaction_id() -> Response:
    """
    Gets the next reaction id for the active workbook

    Returns:
        flask.Response: A JSON response with the new reaction id
    """
    workbook_name = request.json.get("workbook")
    workgroup_name = request.json.get("workgroup")

    if workgroup_name is None or workbook_name is None:
        return jsonify("Bad Request"), status.HTTP_400_BAD_REQUEST

    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )

    new_id = services.reaction.get_next_reaction_id_for_workbook(workbook.id)

    return jsonify(new_id)
