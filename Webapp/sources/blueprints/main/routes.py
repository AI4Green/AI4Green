from flask import (
    Response,
    current_app,
    flash,
    get_template_attribute,
    jsonify,
    redirect,
    render_template,
    request,
    send_file,
    url_for,
    session
)

from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from sources import models, services
from sources.auxiliary import get_notification_number, get_workgroups
from sources.blueprints.auth.forms import LoginForm
from sources.decorators import workbook_member_required
from sources.extensions import db
from flask_jwt_extended import (
    create_access_token,
    create_refresh_token,
)
from werkzeug.urls import url_parse
from datetime import timedelta


from . import main_bp  # imports the blueprint of the main route


# The standard user page is rendered
@main_bp.route("/", methods=["GET", "POST"])
@main_bp.route("/home", methods=["GET", "POST"])
def index() -> Response:
    messages_from_redirects = (
        [request.args.get("message")] if request.args.get("message") else []
    )

    user_confirmed = None
    form = LoginForm()
    user_role = None
    news_items = []

    if request.method == "POST":
        print("Login form submitted from index page")
        print(f"Form data: {form.data}")
        if form.validate_on_submit():
            try:
                # change page_redirect to user.
                user = services.auth.verify_login(form)
                print(f"Verify login result: {user}")
                if user:
                    # After successful login, redirect to the desired page
                    # You can use the same logic as in the /auth/login route
                    user_data = {
                        "id": user.id,
                        "email": user.email,
                        "username": user.username,
                        "role_id": user.role,
                        "Role": user.Role.name,
                        "fullname": user.fullname,
                        "is_verified": user.is_verified,
                        "time_of_creation": user.time_of_creation.isoformat() if user.time_of_creation else None,
                    }

                    access_token = create_access_token(
                        identity=user.id, additional_claims=user_data, expires_delta=timedelta(hours=1))
                    refresh_token = create_refresh_token(
                        identity=user.id, additional_claims=user_data)

                    session['access_token'] = access_token
                    session['refresh_token'] = refresh_token

                    next_page = request.args.get("next")
                    if not next_page or url_parse(next_page).netloc != "":
                        next_page = url_for("main.index")

                    return redirect(next_page)  # Return a redirect.

            except Exception as e:
                print(f"Error during login on index: {e}")
                flash("Something went wrong with login.")
                return redirect(url_for("main.index"))
        else:
            print(f"Form validation failed: {form.errors}")

    if current_user.is_authenticated:
        form = None
        user_role = current_user.Role.name
        user_confirmed = current_user.is_verified

        news_items = (
            db.session.query(models.NewsItem)
            .order_by(models.NewsItem.time.desc())
            .all()
        )
        return render_template(
            "home.html",
            user_role=user_role,
            user_confirmed=user_confirmed,
            news_items=news_items,
            messages_from_redirects=messages_from_redirects,
            form=form,
        )
    else:
        return render_template("landing_page.html", form=form)


@main_bp.route("/load_icons", methods=["GET", "POST"])
def load_icons() -> Response:
    """
    This function renders the icon macro from macros.html for the quick access panel
    """
    selected = request.json.get("input")
    load_type = request.json.get("load_type")
    icon_names = []
    header = None
    bootstrap_icon = ""
    if load_type == "workgroup":
        workbooks = services.workbook.get_workbooks_from_user_group_combination(
            selected
        )
        icon_names = [i.name for i in workbooks]
        header = "Workbooks in " + selected
        load_type = "workbook"
        bootstrap_icon = "bi bi-journal-text"

    elif load_type == "workbook":
        reactions = services.reaction.list_active_in_workbook(
            workbook=selected,
            workgroup=request.json.get("activeWorkgroup"),
            sort_crit="time",
        )
        icon_names = [i.reaction_id for i in reactions[:11]]
        header = "Recent Reactions in " + selected
        load_type = "reaction"
        bootstrap_icon = "bi bi-eyedropper"

    # load macro template with assigned variables
    icon_macro = get_template_attribute("macros.html", "icon_panel")
    return jsonify(icon_macro(icon_names, load_type, header, bootstrap_icon))


@main_bp.route("/get_marvinjs_key", methods=["POST"])
def get_marvinjs_key():
    return jsonify({"marvinjs_key": current_app.config["MARVIN_JS_API_KEY"]})


# Go to the sketcher
@main_bp.route(
    "/sketcher/<workgroup>/<workbook>/<reaction_id>/<tutorial>", methods=["GET", "POST"]
)
@login_required
@workbook_member_required
def sketcher(
    workgroup: str, workbook: str, reaction_id: str, tutorial: str
) -> Response:
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
        reaction={
            "name": "Tutorial Reaction",
            "reaction_id": "TUT-001",
            "reaction_type": "STANDARD",
        },
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
        colours = {
            "Recommended": "#00ff00",
            "Problematic": "#ffff00",
            "Hazardous": "#ff0000",
            "HighlyHazardous": "#8b0000",
            "Recommended_text": "#000000",
            "Problematic_text": "#000000",
            "Hazardous_text": "#000000",
            "HighlyHazardous_text": "#ffffff",
        }
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
