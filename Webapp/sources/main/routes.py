from flask import render_template, request, redirect, url_for, session, flash, send_file, \
    jsonify  # renders html templates
from flask_login import login_required  # protects a view function against anonymous users
from sources.main import main_bp  # imports the blueprint of the main route
from pony.orm import select, desc
from sources import db
from flask_login import current_user
from sources.auxiliary import get_workgroups, get_notification_number, get_workbooks, security_member_workgroup_workbook
from sources import app


# The standard user page is rendered
@main_bp.route('/', methods=['GET', 'POST'])
@main_bp.route('/home', methods=['GET', 'POST'])
@login_required
def index():
    # Set the "role" session variable for a user
    if not session.__contains__("role"):
        session["role"] = select(x.role.name for x in db.User if x.email == current_user.email).first()
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    if request.method == "POST":
        workgroup_selected = request.form['WG-select']
        if workgroup_selected != "-Select Workgroup-":
            return redirect(url_for("workgroup.workgroup", workgroup_selected=workgroup_selected))
    news_items = select(x for x in db.NewsItem).order_by(desc(db.NewsItem.time))[:]
    return render_template('home.html', user_role=session["role"], workgroups=workgroups,
                           notification_number=notification_number, news_items=news_items)


# Go to the sketcher
@main_bp.route('/sketcher/<workgroup>/<workbook>/<reaction_id>/<tutorial>', methods=['GET', 'POST'])
@login_required
def sketcher(workgroup, workbook, reaction_id, tutorial):
    # must be logged in and a member of the workgroup and workbook
    if not security_member_workgroup_workbook(workgroup, workbook):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    reaction = select(rxn for rxn in db.Reaction if rxn.reaction_id == reaction_id and rxn.workbooks.group.name
                      == workgroup).first()
    if reaction.reaction_smiles:
        load_status = "loading"
    else:
        load_status = "loaded"
    return render_template('sketcher_reload.html', reaction=reaction, load_status=load_status,
                           demo="not demo", workgroups=workgroups,
                           notification_number=notification_number, active_workgroup=workgroup,
                           active_workbook=workbook, tutorial=tutorial)


# Go to the sketcher tutorial
@main_bp.route('/sketcher_tutorial/<tutorial>', methods=['GET', 'POST'])
def sketcher_tutorial(tutorial):
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template('sketcher_reload.html', reaction={"name": "Tutorial Reaction", "reaction_id": "TUT-001"},
                           demo="not demo", workgroups=workgroups,
                           notification_number=notification_number, active_workgroup=None,
                           active_workbook=None, tutorial=tutorial)


# Go to demo
@main_bp.route('/demo', methods=['GET', 'POST'])
@login_required
def demo():
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template('demo_sketcher.html', demo="demo", workgroups=workgroups, notification_number=notification_number)


# manage account page
@main_bp.route('/manage_account', methods=['GET', 'POST'])
@login_required
def manage_account():
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template("manage_account.html", workgroups=workgroups, notification_number=notification_number)


# info page
@main_bp.route('/info', methods=['GET', 'POST'])
@login_required
def info():
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template("info.html", workgroups=workgroups, notification_number=notification_number)


# send guide
@main_bp.route('/send_guide', methods=['GET', 'POST'])
@login_required
def send_guide():
    # must be logged in
    return send_file('static/AI4Green_User_Manual.pdf', as_attachment=True)


# send quickstart guide
@main_bp.route('/send_quickstart_guide', methods=['GET', 'POST'])
@login_required
def send_quickstart_guide():
    # must be logged in
    return send_file('static/AI4Green_quick_guide.pdf', as_attachment=True)


# delete reaction
@main_bp.route('/delete_reaction/<reaction_name>/<workgroup>/<workbook>', methods=['GET', 'POST'])
@login_required
def delete_reaction(reaction_name, workgroup, workbook):
    # must be logged in a member of the workgroup and workbook
    if not security_member_workgroup_workbook(workgroup, workbook):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    # find reaction
    reaction = select(x for x in db.Reaction if x.workbooks.name == workbook and x.workbooks.group.name == workgroup and
                      x.creator.user.email == current_user.email and x.name == reaction_name).first()
    # check user is creator of reaction
    if not reaction:
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    # change to inactive
    reaction.status = "inactive"
    return redirect(url_for("workgroup.workgroup", workgroup_selected=workgroup, workbook_selected=workbook))


# marvin js help page
@main_bp.route('/marvin_js_help', methods=['GET', 'POST'])
@login_required
def marvin_js_help():
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template('marvin_js_help.html', workgroups=workgroups, notification_number=notification_number)


# accessibility
@main_bp.route('/accessibility', methods=['GET', 'POST'])
@login_required
def accessibility():
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template('accessibility.html', workgroups=workgroups, notification_number=notification_number)


@main_bp.route('/get_custom_colours', methods=['GET', 'POST'])
@login_required
def get_custom_colours():
    colours = current_user.hazard_colors
    return jsonify({"colours": colours})


@main_bp.route('/change_hazard_colours', methods=['GET', 'POST'])
@login_required
def change_hazard_colours():
    current_user.hazard_colors = {"Recommended": request.form['Recommended'],
                                  "Problematic": request.form['Problematic'],
                                  "Hazardous": request.form['Hazardous'],
                                  "HighlyHazardous": request.form['HighlyHazardous'],
                                  "Recommended_text": request.form['Recommended_text'],
                                  "Problematic_text": request.form['Problematic_text'],
                                  "Hazardous_text": request.form['Hazardous_text'],
                                  "HighlyHazardous_text": request.form['HighlyHazardous_text']}
    if request.form['mode'] == "update":
        return jsonify({"message": "Hazard colours were successfully updated!"})
    else:
        return jsonify({"message": "Hazard colours were reverted to the default!"})


@app.errorhandler(500)
def internal_error(error):
    flash("Something went wrong, please try again later.")
    return redirect(url_for("main.index"))
