from flask import render_template, jsonify, request  # renders html templates
from flask_login import login_required  # protects a view function against anonymous users
from sources.standard_landing_page import standard_landing_page_bp  # imports the blueprint of route
from sources.standard_landing_page.reaction_list import get_scheme_list, get_reaction_list
from sources.auxiliary import get_workgroups, get_notification_number


# Go to the delete profile page
@standard_landing_page_bp.route('/delete_profile', methods=['GET', 'POST'])
@login_required
def delete_profile():
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template('delete_profile.html', workgroups=workgroups, notification_number=notification_number)


# Get reactions
@standard_landing_page_bp.route('/get_reactions', methods=['GET', 'POST'])
@login_required
def get_reactions():
    # must be logged in
    sort_crit = request.form['sort_crit']
    workbook = str(request.form['workbook'])
    workgroup = str(request.form['workgroup'])
    reactions = get_reaction_list(workbook, workgroup, sort_crit)
    reaction_details = render_template('_saved_reactions.html', reactions=reactions, sort_crit=sort_crit)
    return jsonify({'reactionDetails': reaction_details})


@standard_landing_page_bp.route('/get_schemata', methods=['GET', 'POST'])
@login_required
def get_schemata():
    # must be logged in
    workbook = str(request.form['workbook'])
    workgroup = str(request.form['workgroup'])
    size = str(request.form['size'])
    sort_crit = request.form['sort_crit']
    schemes = get_scheme_list(workbook, workgroup, sort_crit, size)
    return {"schemes": schemes, "sort_crit": sort_crit}
