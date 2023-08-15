from flask import (Response, jsonify,  # renders html templates
                   render_template, request)
from flask_login import \
    login_required  # protects a view function against anonymous users

from sources.auxiliary import get_notification_number, get_workgroups

from . import reaction_list_bp  # imports the blueprint of route
from .reaction_list import get_reaction_list, get_scheme_list


# Go to the delete profile page
@reaction_list_bp.route("/delete_profile", methods=["GET", "POST"])
@login_required
def delete_profile() -> Response:
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "delete_profile.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


# Get reactions
@reaction_list_bp.route("/get_reactions", methods=["GET", "POST"])
@login_required
def get_reactions() -> Response:
    # must be logged in
    sort_crit = request.form["sort_crit"]
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
    reaction_id = request.form.get('reaction_id')
    workbook = str(request.form["workbook"])
    workgroup = str(request.form["workgroup"])
    size = str(request.form["size"])
    sort_crit = request.form["sort_crit"]
    schemes = get_scheme_list(workbook, workgroup, sort_crit, size, reaction_id)
    return {"schemes": schemes, "sort_crit": sort_crit}