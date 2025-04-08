"""
This module contains user authentication functions:
login, logout, and register
"""

import pytz

# imports the objects of the login and registration forms
from flask import Response, flash, redirect, render_template, request, session, url_for, Markup
from flask_oidc import OpenIDConnect

# url_parse parses the URL if it is relative or
# absolute to avoid redirection to a malicious site
from flask_login import current_user, login_user, logout_user
from sources import auxiliary, models, services  # imports the user database object

# login_user registers the user as logged in and sets the current_user variable for that user
# logout_user offers users the option to log out of the application
# current_user is a proxy for the current user
from sources.extensions import db
from sqlalchemy import func
from werkzeug.urls import url_parse
from urllib.parse import urlparse
from flask import Response, jsonify
from flask_jwt_extended import (
    create_access_token,
    create_refresh_token,
)
from datetime import timedelta

from . import auth_bp  # imports the blueprint of the

oidc = OpenIDConnect()
from . import auth_bp  # imports the blueprint of the
from .forms import LoginForm, RegistrationForm
from .utils import login_not_allowed



# Login page
@auth_bp.route("/login")
@oidc.require_login
def login():
    """Handles user login via OIDC and retrieves user info."""
    user_info = oidc.user_getinfo(["sub", "email", "name", "preferred_username"])

    if not user_info:
        return jsonify({"error": "Failed to retrieve user info"}), 401

    # Store user info in session (if needed)
    session["user"] = {
        "id": user_info.get("sub"),  # User's unique identifier from OIDC provider
        "email": user_info.get("email"),
        "fullname": user_info.get("name"),
        "username": user_info.get("preferred_username"),
    }

    next_page = request.args.get("next") or url_for("main.index")
    return redirect(next_page)



# Logout option redirecting to the login page
@auth_bp.route("/logout")
def logout():
    """Logs the user out and revokes their OIDC session."""
    oidc.logout()
    session.clear()  # Clear session data
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
    return render_template("privacy_notice.html")


@auth_bp.route("/hazard_disclaimer")
def hazard_disclaimer() -> Response:
    """
    Renders the hazard disclaimer page

    Returns:
        flask.Response The hazard disclaimer page
    """
    return render_template("hazards_disclaimer.html")
