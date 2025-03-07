from datetime import datetime

from flask import Markup, Response, flash, redirect, url_for
from flask_login import current_user
from sources import services
from sources.extensions import db

from . import email_verification_bp


@email_verification_bp.route("/verify_email/<token>")
def verify_email(token: str) -> Response:
    """
    Verifies the email of a user

    Args:
        token: The token to verify the email

    Returns:
        flask.Response: A redirect to the login page with a success or failure message

    """
    user = services.email.verify_encoded_token(token=token, identifier="verify_email")

    if not user:
        flash(
            Markup(
                "Email verification link is invalid or has expired. Please login to resend."
            )
        )
        return redirect(url_for("auth.login"))

    user.is_verified = True
    user.confirmed_on = datetime.now()
    db.session.commit()
    flash("Your account is now verified. Thank you!", "success")

    return redirect(url_for("auth.login"))


@email_verification_bp.route("/email_verification_request", methods=["GET", "POST"])
@email_verification_bp.route(
    "/email_verification_request/<user_id>", methods=["GET", "POST"]
)
def email_verification_request(user_id=None) -> Response:
    """
    Sends an email verification request to a user

    Args:
        user_id:

    Returns:
        flask.Response: A redirect to the login page with a message
    """
    if not user_id:
        user = current_user
    else:
        user = services.user.from_id(int(user_id))

    services.email.send_email_verification(user)
    flash("Thank you! An email has been sent to your inbox.")

    return redirect(url_for("auth.login"))
