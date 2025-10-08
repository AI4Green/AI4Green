"""
This module contains user authentication functions:
login, logout, and register
"""

import uuid

from flask import Response, flash, jsonify, redirect, render_template, request, url_for
from flask_login import current_user, login_user, logout_user
from Webapp.sources.services.user import from_email
from sources import auxiliary, models, services
from sources.extensions import db, oidc

from . import auth_bp
from .forms import LoginForm, RegistrationForm
from .utils import login_not_allowed


@auth_bp.route("/login", methods=["GET", "POST"])
@login_not_allowed
def login() -> Response:
    """
    Handles user authentication and login requests

    Returns:
        flask.Response Either a rendered template landing_page.html
        or a redirect response if the user successfully logs in.
    """
    form = LoginForm()
    if request.method == "POST":
        page_redirect = services.auth.verify_login(form)
    if not current_user.is_authenticated:
        return render_template("auth/login.html", title="Sign In", form=form)
    else:
        return page_redirect


@auth_bp.route("/logout")
def logout() -> Response:
    """
    Logs out the user and redirects to the index page

    Returns:
        flask.Response The landing page
    """
    # Log out of the flask user system
    logout_user()

    # This line has been commented out until SSO feature is fixed
    # Log out of the OIDC system if the user signed in that way.
    # if oidc.user_loggedin:
    #     return oidc.logout(return_to=url_for("main.index"))

    return redirect(url_for("main.index"))


# Registration page


@auth_bp.route("/register", methods=["GET", "POST"])
@login_not_allowed
def register() -> Response:
    """
    Manages user registration

    Returns:
        flask.Response Renders the login template if registration is successful or registration template if not

    """
    form = RegistrationForm()
    if form.validate_on_submit():
        """The form.validate_on_submit returns True when the browser sends the POST
        request as a result of the user pressing the submit button and if all the fields
        passes validation. It returns False when the browser sends the GET request to
        receive the web page with the form or if at least one field fails validation."""
        # Creates a person and user and commits to the database
        p = models.Person()
        # Capitalize unique fields for consistency
        db.session.add(p)
        fullname = auxiliary.sanitise_user_input(form.fullname.data)
        services.user.add(
            form.username.data, form.email.data, fullname, form.password.data, p
        )
        services.email_services.send_email_verification(p.user)

        flash(
            "An email has been sent to your address. Please follow the instructions to verify your account."
        )  # flashes the success message
        return redirect(url_for("auth.login"))  # redirects to the login page
    return render_template(
        "auth/register.html", title="Register", form=form
    )  # renders the registration template


@auth_bp.route("/privacy_notice")
def privacy_notice() -> Response:
    """
    Renders the privacy notice page

    Returns:
        flask.Response The privacy notice page
    """
    return render_template("general/privacy_notice.html")


@auth_bp.route("/hazard_disclaimer")
def hazard_disclaimer() -> Response:
    """
    Renders the hazard disclaimer page

    Returns:
        flask.Response The hazard disclaimer page
    """
    return render_template("general/hazards_disclaimer.html")


# @auth_bp.route("/oidc_login")
# def oidc_login() -> Response:
#     """Redirect the user to the OpenID Connect login page.
#
#     Returns:
#         Response: redirect to the OIDC login page.
#     """
#     return oidc.redirect_to_auth_server(url_for("auth.oidc_callback", _external=True))


# @auth_bp.route("/authorize")
# @oidc.require_login
# def oidc_callback() -> Response:
#     """Callback endpoint for when a user logs in or registers via OIDC.
#
#     When an existing user is logging in, simply redirect them to the main screen.
#     When a new user registers an account via OIDC, add the new user to the AI4Green
#     database.
#
#     Returns:
#         Response: Redirect the to main screen when the user logs in.
#     """
#     # Attempt to find a user in the AI4Green database with an email from the OIDC provider
#     user_info = oidc.user_getinfo(["email", "name"])
#     user = services.user.from_email(user_email=user_info["email"])
#
#     # If the user doesn't exist, add them.
#     if user is None:
#         person = models.Person()
#         db.session.add(person)
#         services.user.add(
#             username=user_info["name"],
#             email=user_info["email"],
#             fullname=user_info["name"],
#             # generate UUID for mandatory field
#             password_data=str(uuid.uuid4()),
#             person=person,
#         )
#         # send verification email
#         services.email_services.send_email_verification(person.user)
#         # get the new user
#         user = services.user.from_email(user_email=user_info["email"])
#
#     # OIDC and regular login are different. Use `login_user` to hook into
#     # the regular login/out system
#     login_user(user=user)
#
#     return redirect(url_for("main.index"))


@auth_bp.route("/retrosynthesis_key", methods=["POST"])
def get_retrosynthesis_key() -> Response:
    """Endpoint for getting a user's retrosynthesis key
    based on their email and password.

    Returns:
        Response: The user's retrosynthesis key if they have one.
    """
    data = request.get_json()

    if not "email" in data or "password" not in data:
        return jsonify({"error": "Username and password are required."}), 400

    user = from_email(data["email"])
    if not user:
        return jsonify({"error": "Invalid username or password."}), 401

    if not user.check_password(data["password"]):
        return jsonify({"error": "Invalid username or password."}), 401

    if user.retrosynthesis_access_key:
        return jsonify({"access_key": user.retrosynthesis_access_key.key}), 200
    else:
        return jsonify({"message": "User has no access key."}), 200
