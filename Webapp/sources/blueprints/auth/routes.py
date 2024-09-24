"""
This module contains user authentication functions:
login, logout, and register
"""
from datetime import datetime

import pytz

# imports the objects of the login and registration forms
from flask import Response, flash, redirect, render_template, request, session, url_for, Markup

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

from . import auth_bp  # imports the blueprint of the
from .forms import LoginForm, RegistrationForm
from .utils import login_not_allowed

# render_template renders html templates
# redirect instructs the client web browser to automatically
# navigate to a different page, given as an argument
# url_for generates URLs using its internal mapping of URLs to view functions
# flash stores a message for the user to show it when called from the templates
# request parses incoming request data and gives access to it


# Login page
@auth_bp.route("/login", methods=["GET", "POST"])
@login_not_allowed
def login() -> Response:  # the login view function
    # anyone may view

    form = LoginForm()  # instantiates an object of LoginForm
    if request.method == "POST":
        page_redirect = services.auth.verify_login(form)
    if not current_user.is_authenticated:
        return render_template(
            "auth/login.html", title="Sign In", form=form
        )  # renders the login template
    else:
        return page_redirect


# Logout option redirecting to the login page
@auth_bp.route("/logout")
def logout() -> Response:  # the logout view function
    # anyone may view
    logout_user()
    return redirect(url_for("main.index"))


# Registration page


@auth_bp.route("/register", methods=["GET", "POST"])
@login_not_allowed
def register() -> Response:  # the view function that handles user registrations
    # anyone may view
    form = RegistrationForm()  # instantiates an object of RegistrationForm
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
    # anyone may view
    return render_template("privacy_notice.html")


@auth_bp.route("/hazard_disclaimer")
def hazard_disclaimer() -> Response:
    # anyone may view
    return render_template("hazards_disclaimer.html")
