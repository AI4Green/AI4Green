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
from flask import Response, jsonify
from flask_jwt_extended import (
    create_access_token,
    create_refresh_token,
    get_jwt_identity,
    jwt_required,
    unset_jwt_cookies,
)
from datetime import timedelta


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
def login() -> Response:
    form = LoginForm()
    if request.method == "GET":
        return render_template("auth/login.html", title="Sign In", form=form)  # Serve login page
    
    if form.validate_on_submit():
        user = services.auth.verify_login(form)
        print(f"User returned from verify_login: {user}")
        if user:
            # Prepare user data to be added to the token payload
            user_data = {
                "id": user.id,
                "email": user.email,
                "username": user.username,
                "role_id": user.role,
                "Role": user.Role.name,
                "fullname": user.fullname,  # Add other user data as needed
                "is_verified": user.is_verified,
                "time_of_creation": user.time_of_creation.isoformat() if user.time_of_creation else None,
            }

            # Create access and refresh tokens with user data in the payload
            access_token = create_access_token(identity=user.id, additional_claims=user_data, expires_delta=timedelta(hours=1))
            refresh_token = create_refresh_token(identity=user.id, additional_claims=user_data)

            print(f"Access token: {access_token}")
            print(f"Refresh token: {refresh_token}")

            # You can choose to store the token in the session or cookies
            session['access_token'] = access_token
            session['refresh_token'] = refresh_token

            # Determine the next page
            next_page = request.args.get("next")
            if not next_page or url_parse(next_page).netloc != "":
                next_page = url_for("main.index")  # Default redirect URL

            # Redirect to the next page after login
            return redirect(next_page)  # This will redirect to the main index

    return jsonify({"error": "Invalid credentials"}), 401  # Return error message if login fails



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
