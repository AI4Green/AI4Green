from flask import Response, flash, redirect, render_template, request, url_for
from flask_login import current_user, login_required
from sources import models, services
from sources.auxiliary import get_notification_number, get_workgroups
from sources.extensions import db

from . import update_email_bp
from .forms import UpdateEmailForm


@update_email_bp.route("/update_email_password", methods=["GET", "POST"])
@login_required
@update_email_bp.doc(security="sessionAuth")
def update_email_password() -> Response:
    """
    Update email and password page

    Returns:
        flask.Response: renders the update email and password page
    """
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
            services.email_services.send_password_reset(user)
            flash(
                "Please check your email for instructions on how to change your password"
            )
            return redirect(url_for("main.index"))
    return render_template(
        "account_management/update_email_password.html",
        email=current_email,
    )


@update_email_bp.route("/update_email", methods=["GET", "POST"])
@login_required
@update_email_bp.doc(security="sessionAuth")
def update_email() -> Response:
    """
    Update email and sends verification email to the new email address.

    Returns:
        flask.Response: renders the update email page.
    """
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
        user.is_verified = False
        user.update()
        services.email_services.send_email_verification(user)
        # alert user of success and redirect to home
        flash(
            "Your email address has been successfully updated! Please check your inbox for instructions on how to verify this new address!"
        )
        return redirect(url_for("main.index"))
    return render_template(
        "account_management/update_email.html",
        form=form,
    )
