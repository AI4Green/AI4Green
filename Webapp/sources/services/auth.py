from typing import Dict

from flask import flash, redirect, request, url_for, Markup, session
from flask_login import current_user, login_user
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook
from sources.extensions import db
from werkzeug.urls import url_parse
from urllib.parse import urlparse
from sqlalchemy import func
from datetime import datetime


def verify_login(form) -> redirect:
    """
    Logic to validate login request.
    Args:
        form: Form with login credentials

    Returns:
        redirect to home page or previous page that required login
    """
    if form.validate_on_submit():
        """The form.validate_on_submit returns True when the browser sends the POST
        request as a result of the user pressing the submit button and if all the fields
        passes validation. It returns False when the browser sends the GET request to
        receive the web page with the form or if at least one field fails validation."""
        input_data = form.username.data.lower()
        user = (
            db.session.query(models.User)
            .filter(
                (func.lower(models.User.username) == input_data) |
                (func.lower(models.User.email) == input_data)
            )
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
            return redirect(url_for("main.index"))

        # if user was added after 22/04/2024 their email needs to be verified before login
        if user.time_of_creation and user.time_of_creation > datetime(2024, 4, 22) and not user.is_verified:
            verification_url = "/email_verification_request/" + str(user.id)
            flash(
                Markup(
                    f"Please verify your email address before logging in. Didn't receive an email? <a href='{verification_url}' class='alert-link'>Click here</a> to resend."
                )
            )
            return redirect(url_for("main.index"))

        login_user(user, remember=form.remember_me.data)
        role = (
            db.session.query(models.Role)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .first()
        )
        user.most_recent_login_location = services.utils.get_location()
        db.session.commit()
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
            return redirect(url_for("main.index"))

        else:
            """If the login URL includes a next argument that is set
            to a relative path (a URL without the domain portion),
            then the user is redirected to that URL."""
            next_page = next_page.replace('\\', '')
            if not urlparse(next_page).netloc:
                return redirect(next_page)


def reaction_files(
    permission_level: str,
    request_source: str = "user",
    file_object_for_deletion: models.ReactionDataFile = None,
):
    """
    Authenticates user to either view or edit the reaction.
    Args:
        permission_level: Takes value of 'edit' or 'view_only'
        request_source: The origin of the request, either 'user' or 'server'
        file_object_for_deletion: If deleting a file object, we pass that object as an argument.
    """
    if permission_level == "view_only":
        view_files()
    elif permission_level == "edit":
        reaction = services.reaction.get_current_from_request()
        edit_reaction(
            reaction,
            request_source,
            file_attachment=True,
            file_object_for_deletion=file_object_for_deletion,
        )


def view_files():
    """Authenticates user as a workbook member or aborts. Gets the workgroup_name, workbook_name, and workbook."""
    if request.method == "GET":
        workgroup_name = request.args.get("workgroup")
        workbook_name = request.args.get("workbook")
    else:
        workgroup_name = request.form["workgroup"]
        workbook_name = request.form["workbook"]
    # validate user belongs to the workbook
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)


def edit_reaction(
    reaction: models.Reaction,
    request_source: str = "user",
    file_attachment: bool = False,
    file_object_for_deletion: models.ReactionDataFile = None,
):
    """
    Validates the active user has permission to edit the reaction. Protects against user edited HTML.
    Aborts process with a 401 error if validation is failed

    Args:
        reaction: the active reaction
        request_source: the origin of request. Either 'user' or 'server'
        file_attachment: Whether file attachments are being edited.
        file_object_for_deletion: If a file is being deleted, we check this file can be deleted
    """
    # validate user is in workbook
    workbook_persons = reaction.workbook.users
    workbook_users = [x.user for x in workbook_persons]
    if current_user not in workbook_users:
        abort(401)
    # validate the user is the creator
    if reaction.creator_person.user.email != current_user.email:
        abort(401)
    # validate the user is not deleting autogenerated files.
    if (
        file_attachment
        and request.method == "DELETE"
        and request_source == "user"
        and file_object_for_deletion.autogenerated is True
    ):
        abort(401)
    if file_attachment:
        return
    # validate the reaction is not locked, unless it is a file attachment being edited.
    if reaction.complete == "complete":
        abort(401)
