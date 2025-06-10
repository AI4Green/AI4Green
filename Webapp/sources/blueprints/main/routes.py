from datetime import datetime
from typing import Literal

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
)
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from sources import models, services
from sources.auxiliary import (
    get_notification_number,
    get_workgroups,
    security_pi_workgroup,
)
from sources.blueprints.auth.forms import LoginForm
from sources.decorators import principal_investigator_required, workbook_member_required
from sources.extensions import db

from . import main_bp


@main_bp.route("/", methods=["GET", "POST"])
@main_bp.route("/home", methods=["GET", "POST"])
def index() -> Response:
    """
    The home page is rendered if the user is authenticated, otherwise the landing page is rendered.

    Returns:
        flask.Response:  Either a rendered template (home.html or landing_page.html)
        or a redirect response from the auth service.
    """
    # used to display flash messages after a redirect following a fetch request.
    messages_from_redirects = (
        [request.args.get("message")] if request.args.get("message") else []
    )

    form = LoginForm()
    privacy_policy_accepted = True
    privacy_policy_date = services.utils.get_privacy_policy_date()

    if request.method == "POST":
        # return redirects from login verification to prevent form resubmission
        page_redirect = services.auth.verify_login(form)
        return page_redirect

    # update jinja variables if user is logged in to populate the homepage
    if current_user.is_authenticated:
        form = None
        user_role = current_user.Role.name
        user_confirmed = current_user.is_verified
        user_privacy_policy = current_user.privacy_policy_accepted_on

        if user_privacy_policy is None or user_privacy_policy < privacy_policy_date:
            privacy_policy_accepted = False

        news_items = (
            db.session.query(models.NewsItem)
            .order_by(models.NewsItem.time.desc())
            .all()
        )
        return render_template(
            "general/home.html",
            user_role=user_role,
            user_confirmed=user_confirmed,
            news_items=news_items,
            messages_from_redirects=messages_from_redirects,
            form=form,
            privacy_policy_accepted=privacy_policy_accepted,
        )
    # user is not authenticated, send to landing page.
    else:
        return render_template("general/landing_page.html", form=form)


@main_bp.route("/load_icons", methods=["GET", "POST"])
def load_icons() -> Response:
    """
    Renders the icon macro from macros.html for the quick access panel.

    Returns:
        flask.Response: A JSON response containing the rendered icon macro.
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
    icon_macro = get_template_attribute("macros/icons.html", "icon_panel")
    return jsonify(icon_macro(icon_names, load_type, header, bootstrap_icon))


@main_bp.route("/get_marvinjs_key", methods=["POST"])
@main_bp.doc(hide=True)
def get_marvinjs_key():
    """
    Retrieves the marvin JS API key from the configuration file.
    Returns:
        flask.Response: A JSON response containing the marvin JS API key.
    """
    return jsonify({"marvinjs_key": current_app.config["MARVIN_JS_API_KEY"]})


@main_bp.route(
    "/sketcher/<workgroup>/<workbook>/<reaction_id>/<tutorial>", methods=["GET", "POST"]
)
@login_required
@workbook_member_required
@main_bp.doc(
    security="sessionAuth",
    description="Requires login and membership in the specified workbook.",
)
def sketcher(
    workgroup: str, workbook: str, reaction_id: str, tutorial: Literal["yes", "no"]
) -> Response:
    """
    Renders the sketcher page with the given reaction ID.

    Args:
        workgroup: the workgroup the reaction belongs to
        workbook: the workbook the reaction belongs to
        reaction_id: the reaction ID of the reaction to be loaded
        tutorial: whether the user is in tutorial mode

    Returns:
        flask.Response The rendered sketcher page.
    """
    workbook_object = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup, workbook
    )
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook_object.id
    )
    addenda = services.reaction.get_addenda(reaction)

    if reaction.reaction_smiles:
        load_status = "loading"
    else:
        load_status = "loaded"
    return render_template(
        "reactions/sketcher_reload.html",
        reaction=reaction,
        load_status=load_status,
        demo="not demo",
        active_workgroup=workgroup,
        active_workbook=workbook,
        tutorial=tutorial,
        addenda=addenda,
        review=False,
    )


@main_bp.route("/sketcher_tutorial/<tutorial>", methods=["GET", "POST"])
def sketcher_tutorial(tutorial: str) -> Response:
    """
    Renders the tutorial sketcher page

    Args:
        tutorial: whether the user is in tutorial mode

    Returns:
        flask.Response: The rendered tutorial sketcher page.
    """
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "reactions/sketcher_reload.html",
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


@main_bp.route("/demo", methods=["GET", "POST"])
def demo() -> Response:
    """
    Renders the demo sketcher page.

    Returns:
        flask.Response: The rendered demo sketcher page.
    """
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "reactions/demo_sketcher.html",
        demo="demo",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@main_bp.route("/search", methods=["GET", "POST"])
@login_required
@main_bp.doc(security="sessionAuth")
def search() -> Response:
    """
    Renders the search page.

    Returns:
        flask.Response: The rendered search page.
    """
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "reactions/search.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@main_bp.route("/accept_privacy_policy", methods=["PATCH"])
def accept_privacy_policy() -> Response:
    """
    Accepts the privacy policy and updates the property in the database for the current user.

    Returns:
        flask.Response: 204 successful with no additional content
    """
    user = services.user.from_email(current_user.email)
    user.privacy_policy_accepted_on = datetime.now()
    db.session.commit()
    return Response(status=204)


# manage account page
@main_bp.route("/manage_account", methods=["GET", "POST"])
@login_required
@main_bp.doc(security="sessionAuth")
def manage_account() -> Response:
    """
    Renders the manage account page.

    Returns:
        flask.Response: The rendered manage account page.
    """
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "account_management/manage_account.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


# info page
@main_bp.route("/info", methods=["GET", "POST"])
def info() -> Response:
    """
    Renders the info (help) page.

    Returns:
        flask.Response: The rendered info page.
    """
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "general/info.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


# about page
@main_bp.route("/about", methods=["GET", "POST"])
def about() -> Response:
    """
    Renders the about page

    Returns:
        flask.Response: The rendered about page.
    """
    workgroups = []
    notification_number = 0
    if current_user.is_authenticated:
        workgroups = get_workgroups()
        notification_number = get_notification_number()
    return render_template(
        "general/about.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


# send guide
@main_bp.route("/send_guide", methods=["GET", "POST"])
def send_guide() -> Response:
    """
    Sends the user manual as a file download.

    Returns:
        flask.Response: The user manual file.
    """
    return send_file("static/AI4Green_User_Manual.pdf", as_attachment=True)


# send quickstart guide
@main_bp.route("/send_quickstart_guide", methods=["GET", "POST"])
def send_quickstart_guide() -> Response:
    """
    Sends the quick start guide as a file download.

    Returns:
        flask.Response: The quick start guide file.
    """
    return send_file("static/AI4Green_quick_guide.pdf", as_attachment=True)


# marvin js help page
@main_bp.route("/marvin_js_help", methods=["GET", "POST"])
def marvin_js_help() -> Response:
    """
    Renders the Marvin JS help page.

    Returns:
        flask.Response: The rendered Marvin JS help page.
    """
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "general/marvin_js_help.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@main_bp.route("/accessibility", methods=["GET", "POST"])
@login_required
@main_bp.doc(security="sessionAuth")
def accessibility() -> Response:
    """
    Renders the accessibility page where the sustainability colour-coding can be changed.

    Returns:
        flask.Response: The rendered accessibility page.
    """
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "account_management/accessibility.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@main_bp.route("/get_custom_colours", methods=["GET", "POST"])
def get_custom_colours() -> Response:
    """
    Retrieves the custom colours for the sustainability colour-coding. Default if user is not logged in.

    Returns:
        flask.Response: A JSON response containing the custom colours.
    """
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
@main_bp.doc(security="sessionAuth")
def change_hazard_colours() -> Response:
    """
    Changes the sustainability colour-coding for the user.

    Returns:
        flask.Response: A JSON response indicating the success of the operation.
    """
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
