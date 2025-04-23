"""
This module contains user authentication functions:
login, logout, and register
"""

from flask import Response, flash, redirect, render_template, request, url_for
from flask_login import current_user, logout_user
from sources import auxiliary, models, services
from sources.extensions import db

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
    logout_user()
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
        services.email.send_email_verification(p.user)

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
