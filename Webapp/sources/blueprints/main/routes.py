from flask import (Response, flash, redirect, render_template, request,
                   send_file, url_for, jsonify, current_app)
from flask_login import (  # protects a view function against anonymous users
    current_user, login_required)

from sources import models
from sources.auxiliary import (get_notification_number, get_workbooks,
                               get_workgroups,
                               security_member_workgroup_workbook)
from sources.extensions import db

from . import main_bp  # imports the blueprint of the main route


# The standard user page is rendered
@main_bp.route("/", methods=["GET", "POST"])
@main_bp.route("/home", methods=["GET", "POST"])
def index() -> Response:
    # gets jinja variables if user is logged in to populate the homepage, else the non-logged in user homepage is loaded
    if current_user.is_authenticated:
        # get the users role to determine if they can access the admin dashboard button or not
        user = (
            db.session.query(models.User)
            .filter(models.User.email == current_user.email)
            .first()
        )
        user_role = user.Role.name
        workgroups = get_workgroups()
        notification_number = get_notification_number()
        if request.method == "POST":
            workgroup_selected = request.form["WG-select"]
            if workgroup_selected != "-Select Workgroup-":
                return redirect(
                    url_for("workgroup.workgroup", workgroup_selected=workgroup_selected)
                )
        news_items = (
            db.session.query(models.NewsItem).order_by(models.NewsItem.time.desc()).all()
        )
    else:
        # user not logged in
        user_role = None
        workgroups = []
        notification_number = 0
        news_items = []
    return render_template(
        "home.html",
        user_role=user_role,
        workgroups=workgroups,
        notification_number=notification_number,
        news_items=news_items,
    )


@main_bp.route(
    "/get_marvinjs_key", methods=["POST"]
)
def get_marvinjs_key():
    return jsonify({'marvinjs_key': current_app.config["MARVIN_JS_API_KEY"]})


# Go to the sketcher
@main_bp.route(
    "/sketcher/<workgroup>/<workbook>/<reaction_id>/<tutorial>", methods=["GET", "POST"]
)
@login_required
def sketcher(
    workgroup: str, workbook: str, reaction_id: str, tutorial: str
) -> Response:
    # validates the user is a member of the workgroup and workbook by checking the database
    if not security_member_workgroup_workbook(workgroup, workbook):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    workbook_object = (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )
    reaction = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .filter(models.WorkBook.id == workbook_object.id)
        .first()
    )
    addenda = (
        db.session.query(models.ReactionNote)
        .join(models.Reaction)
        .filter(models.Reaction.id == reaction.id)
        .all()
    )
    file_attachments = reaction.file_attachments
    if reaction.file_attachments:
        file_attachments = reaction.file_attachments
    if reaction.reaction_smiles:
        load_status = "loading"
    else:
        load_status = "loaded"
    return render_template(
        "sketcher_reload.html",
        reaction=reaction,
        load_status=load_status,
        demo="not demo",
        workgroups=workgroups,
        notification_number=notification_number,
        active_workgroup=workgroup,
        active_workbook=workbook,
        tutorial=tutorial,
        addenda=addenda,
        file_attachments=file_attachments,
    )


# Go to the sketcher tutorial
@main_bp.route("/sketcher_tutorial/<tutorial>", methods=["GET", "POST"])
def sketcher_tutorial(tutorial: str) -> Response:
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "sketcher_reload.html",
        reaction={"name": "Tutorial Reaction", "reaction_id": "TUT-001"},
        demo="not demo",
        workgroups=workgroups,
        notification_number=notification_number,
        active_workgroup=None,
        active_workbook=None,
        tutorial=tutorial,
    )


# Go to demo
@main_bp.route("/demo", methods=["GET", "POST"])
def demo() -> Response:
    # must be logged in
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "demo_sketcher.html",
        demo="demo",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@main_bp.route("/search", methods=["GET", "POST"])
@login_required
def search() -> Response:
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "search.html", workgroups=workgroups, notification_number=notification_number
    )


# manage account page
@main_bp.route("/manage_account", methods=["GET", "POST"])
@login_required
def manage_account() -> Response:
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "manage_account.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


# info page
@main_bp.route("/info", methods=["GET", "POST"])
def info() -> Response:
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "info.html", workgroups=workgroups, notification_number=notification_number
    )


# about page
@main_bp.route("/about", methods=["GET", "POST"])
def about() -> Response:
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "about.html", workgroups=workgroups, notification_number=notification_number
    )


# send guide
@main_bp.route("/send_guide", methods=["GET", "POST"])
def send_guide() -> Response:
    # must be logged in
    return send_file("static/AI4Green_User_Manual.pdf", as_attachment=True)


# send quickstart guide
@main_bp.route("/send_quickstart_guide", methods=["GET", "POST"])
def send_quickstart_guide() -> Response:
    # must be logged in
    return send_file("static/AI4Green_quick_guide.pdf", as_attachment=True)


# marvin js help page
@main_bp.route("/marvin_js_help", methods=["GET", "POST"])
def marvin_js_help() -> Response:
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "marvin_js_help.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


# accessibility
@main_bp.route("/accessibility", methods=["GET", "POST"])
@login_required
def accessibility() -> Response:
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "accessibility.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@main_bp.route("/get_custom_colours", methods=["GET", "POST"])
def get_custom_colours() -> Response:
    if current_user.is_authenticated:
        colours = current_user.hazard_colors
    else:
        # use default colours if user is not logged in
        colours = {'Recommended': '#00ff00',
                   'Problematic': '#ffff00',
                   'Hazardous': '#ff0000',
                   'HighlyHazardous': '#8b0000',
                   'Recommended_text': '#000000',
                   'Problematic_text': '#000000',
                   'Hazardous_text': '#000000',
                   'HighlyHazardous_text': '#ffffff'}
    return jsonify({"colours": colours})


@main_bp.route("/change_hazard_colours", methods=["GET", "POST"])
@login_required
def change_hazard_colours() -> Response:
    current_user.hazard_colors = {
        "Recommended": request.form["Recommended"],
        "Problematic": request.form["Problematic"],
        "Hazardous": request.form["Hazardous"],
        "HighlyHazardous": request.form["HighlyHazardous"],
        "Recommended_text": request.form["Recommended_text"],
        "Problematic_text": request.form["Problematic_text"],
        "Hazardous_text": request.form["Hazardous_text"],
        "HighlyHazardous_text": request.form["HighlyHazardous_text"],
    }
    current_user.update()
    if request.form["mode"] == "update":
        return jsonify({"message": "Hazard colours were successfully updated!"})
    else:
        return jsonify({"message": "Hazard colours were reverted to the default!"})


@main_bp.errorhandler(500)
def internal_error(error):
    flash("Something went wrong, please try again later.")
    return redirect(url_for("main.index"))
