from flask import redirect  # renders html templates
from flask import Response, flash, render_template, request, url_for
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from sources import models, services
from sources.auxiliary import get_notification_number, get_workgroups
from sources.extensions import db

from . import update_email_bp  # imports the blueprint of the dummy route
from .forms import UpdateEmailForm


# Go to the update email/password page
@update_email_bp.route("/update_email_password", methods=["GET", "POST"])
@login_required
def update_email_password() -> Response:
    # must be logged in
    current_email = current_user.email
    user = (
        db.session.query(models.User)
        .filter(models.User.username == current_user.username)
        .first()
    )
    if request.method == "POST":
        if request.form.get("ChangeEmail") == "Value1":
            return redirect(url_for("update_email.update_email"))
        elif request.form.get("UpdatePassword") == "Value2":
            services.email.send_password_reset(user)
            flash(
                "Please check your email for instructions on how to change your password"
            )
            return redirect(url_for("main.index"))
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "update_email_password.html",
        email=current_email,
        workgroups=workgroups,
        notification_number=notification_number,
    )


# Go to the update email page
@update_email_bp.route("/update_email", methods=["GET", "POST"])
@login_required
def update_email() -> Response:
    # must be logged in
    form = UpdateEmailForm()
    if form.validate_on_submit():
        # if authentication is successful
        user = (
            db.session.query(models.User)
            .filter(models.User.username == current_user.username)
            .first()
        )
        if user is None or not user.check_password(form.old_password.data):
            flash("Invalid password")
            return redirect(url_for("update_email.update_email"))
        # check email is not already registered
        emails = (
            db.session.query(models.User)
            .filter(models.User.email == form.email.data)
            .first()
        )
        if emails:
            flash("This email is already registered to an account")
            return redirect(url_for("update_email.update_email"))
        # change user's email
        user.email = form.email.data
        user.update()
        # alert user of success and redirect to home
        flash("Your email has been successfully updated!")
        return redirect(url_for("main.index"))
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "update_email.html",
        form=form,
        workgroups=workgroups,
        notification_number=notification_number,
    )
