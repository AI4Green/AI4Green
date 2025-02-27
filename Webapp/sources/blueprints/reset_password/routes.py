from flask import redirect  # renders html templates
from flask import Response, flash, render_template, url_for
from flask_login import current_user, logout_user
from sources import models, services
from sources.extensions import db  # imports the user database models

from . import reset_password_bp
from .forms import ResetPasswordForm, ResetPasswordRequestForm


@reset_password_bp.route("/reset_password_request", methods=["GET", "POST"])
def reset_password_request() -> Response:
    """
    Reset password request route

    Returns:
        flask.Response: renders the reset password request template and validates the form and sends reset  email
        or redirects to log in if the user is already logged in

    """
    # anyone can access this
    # if the user is already logged in redirect to home
    if current_user.is_authenticated:
        return redirect(url_for("main.index"))
    # create the form
    form = ResetPasswordRequestForm()
    # when it is validly submitted
    if form.validate_on_submit():
        # find the corresponding email
        user = (
            db.session.query(models.User)
            .filter(models.User.email == form.email.data)
            .first()
        )
        # if found send email and redirect to log in page
        if user:
            services.email.send_password_reset(user)
            flash(
                "Please check your email for instructions on how to change your password"
            )
            return redirect(url_for("auth.login"))
        # else tell the user that email not registered
        else:
            flash("We could not find an account with this email")
    return render_template(
        "reset_password_request.html", title="Reset Password", form=form
    )


@reset_password_bp.route("/reset_password/<token>", methods=["GET", "POST"])
def reset_password(token: str) -> Response:
    """
    Reset password route from email link

    Args:
        token: token from email link

    Returns:
        flask.Response: renders the reset password template and validates the form and updates the password

    """
    # token must be correct
    # log out user if they are logged in
    logout_user()
    # if the user is already logged in redirect to home
    if current_user.is_authenticated:
        return redirect(url_for("main.index"))
    # verify link token is correct
    user = services.email.verify_encoded_token(token=token, identifier="reset_password")
    # if link has expired send back to log in
    if not user:
        flash("Password reset link expired")
        return redirect(url_for("auth.login"))
    # create reset form
    form = ResetPasswordForm()
    # when validly submitted
    if form.validate_on_submit():
        # update password and send back to log in
        user.password_hash = user.set_password(form.password.data)
        user.update()
        flash("Your password has been reset!")
        return redirect(url_for("auth.login"))
    return render_template("reset_password.html", form=form)
