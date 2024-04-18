from flask import redirect, get_flashed_messages  # renders html templates
from flask import Response, flash, render_template, url_for
from flask_login import current_user, logout_user, login_required
from datetime import datetime
from sources import models, services
from sources.extensions import db  # imports the user database models

from . import email_verification_bp  # imports the blueprint of the dummy route


@email_verification_bp.route("/verify_email/<token>")
def verify_email(token: str) -> Response:
    if current_user.is_confirmed:
        flash("Your account has already been confirmed!")
        return redirect(url_for("auth.login"))
    user = services.email.verify_reset_password_token(token)
    if not user:
        flash("Email verification link is invalid or has expired.")
        return redirect(url_for("auth.login"))

    user.is_confirmed = True
    user.confirmed_on = datetime.now()
    db.session.commit()
    flash("Your account is now verified. Thank you!", "success")

    return redirect(url_for("auth.login"))


@email_verification_bp.route("/email_verification_request", methods=["GET", "POST"])
def email_verification_request() -> Response:
    services.email.send_email_verification(current_user)
    flash('Thank you! An email has been sent to your inbox.')

    return redirect(url_for("auth.login"))
