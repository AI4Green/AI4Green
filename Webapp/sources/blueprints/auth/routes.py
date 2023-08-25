"""
This module contains user authentication functions:
login, logout, and register
"""
# imports the objects of the login and registration forms
from flask import Response, flash, redirect, render_template, request, session, url_for

# url_parse parses the URL if it is relative or
# absolute to avoid redirection to a malicious site
from flask_login import current_user, login_user, logout_user
from sources import auxiliary  # imports the user database object
from sources import models

# login_user registers the user as logged in and sets the current_user variable for that user
# logout_user offers users the option to log out of the application
# current_user is a proxy for the current user
from sources.extensions import db
from sqlalchemy import func

# render_template renders html templates
# redirect instructs the client web browser to automatically
# navigate to a different page, given as an argument
# url_for generates URLs using its internal mapping of URLs to view functions
# flash stores a message for the user to show it when called from the templates
# request parses incoming request data and gives access to it

from werkzeug.urls import url_parse
from datetime import datetime
import pytz
from . import auth_bp  # imports the blueprint of the
from .forms import LoginForm, RegistrationForm
from .utils import login_not_allowed


# Login page
@auth_bp.route("/login", methods=["GET", "POST"])
@login_not_allowed
def login() -> Response:  # the login view function
    # anyone may view

    form = LoginForm()  # instantiates an object of LoginForm
    if form.validate_on_submit():
        """The form.validate_on_submit returns True when the browser sends the POST
        request as a result of the user pressing the submit button and if all the fields
        passes validation. It returns False when the browser sends the GET request to
        receive the web page with the form or if at least one field fails validation."""
        user = (
            db.session.query(models.User)
            .filter(func.lower(models.User.username) == form.username.data.lower())
            .first()
        )
        """The select function will search through all of the User entities in the
        database and will return a query that only includes the objects that have
        a matching username. Since there is only one or zero results, the query is
        completed by calling first(), which will return the user object if it exists,
        or None if it does not."""
        if user is None or not user.check_password(form.password.data):
            """If it got a match for the username that was provided, it can next check
            if the password came with the form is valid. This is done by invoking the
            check_password() method defined in models.py. This will take the password
            hash stored with the user and determine if the password entered in the
            form matches the hash or not. In either of two possible error conditions -
            the invalid username or the incorrect password - the error message is
            flashed, and the user is redirected back to the login prompt to try again.
            """
            flash("Invalid username or password")
            return redirect(url_for("auth.login"))
        login_user(user, remember=form.remember_me.data)
        role = (
            db.session.query(models.Role)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .first()
        )
        session["role"] = role.name
        """If the username and password are both correct, then the login_user() function
        from Flask-Login is called. This function will register the user as logged in,
        which means that any future pages the user navigates to will have the current_user
        variable set to that user."""
        next_page = request.args.get("next")
        """The next query string argument is set to the original URL,
        so the application can use that to redirect back after login."""
        if not next_page or url_parse(next_page).netloc != "":
            """If the login URL does not have a next argument or the
            next argument is set to a full URL that includes a domain
            name, then the user is redirected to the index page."""
            next_page = url_for("main.index")
            """If the login URL includes a next argument that is set
            to a relative path (a URL without the domain portion),
            then the user is redirected to that URL."""
        return redirect(next_page)
    return render_template(
        "auth/login.html", title="Sign In", form=form
    )  # renders the login template


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
        fullname = auxiliary.sanitise_user_input(form.fullname.data)
        # Capitalize unique fields for consistency
        db.session.add(p)
        user = models.User(
            username=form.username.data,
            email=form.email.data,
            fullname=fullname,
            Person=p,
            password_hash=models.User.set_password(form.password.data),
        )
        db.session.add(user)
        db.session.commit()

        flash(
            "Congratulations, you are now a registered user!"
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
