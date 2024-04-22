from time import time
from typing import Tuple, Dict

import jwt
from flask import current_app, render_template
from sources import models
from sources.extensions import mail


def get_encoded_token(time_limit: int, arg_dict: Dict) -> str:
    """
    Get token with expiry time.

    Args:
        time_limit: number of seconds before token expires. Use 600 for password resets and 7200 for email verification
        arg_dict: arguments to encode. Should contain either "password_reset" or "email_verification" key with user.id value

    Returns:
        A token string.
    """
    arg_dict["exp"] = time() + time_limit
    return jwt.encode(
        arg_dict,
        current_app.config["SECRET_KEY"],
        algorithm="HS256",
    )


def send_email_verification(user: models.User) -> None:
    """
    Send email verification to defined user.

    Args:
        user: User to send to.
    """
    token = get_encoded_token(time_limit=7200, arg_dict={"verify_email": user.id})
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Green Email Verification",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/email_verification.txt", user=user, token=token, protocol=protocol
        ),
        html_body=render_template(
            "email/email_verification.html", user=user, token=token, protocol=protocol
        ),
    )


def send_password_reset(user: models.User) -> None:
    """
    Send password reset email to a given user.

    Args:
        user: User to send to.
    """
    token = get_encoded_token(time_limit=600, arg_dict={"reset_password": user.id})
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Green Reset Your Password",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/reset_password_text.txt", user=user, token=token, protocol=protocol
        ),
        html_body=render_template(
            "email/reset_password_text.html", user=user, token=token, protocol=protocol
        ),
    )


def send_notification(person: models.Person):
    """
    Send notifications email to a user.

    Args:
        person: Person to send to.

    """
    protocol = get_protocol_type()
    mail.send_email(
        "You have a new AI4Green notification",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[person.user.email],
        text_body=render_template(
            "email/notification_text.txt", user=person.user, protocol=protocol
        ),
        html_body=render_template(
            "email/notification_text.html", user=person.user, protocol=protocol
        ),
    )


def send_password_reset_test(user: models.User) -> Tuple[str, str]:
    """
    Send password reset email in testing.

    Args:
        user: User to send to.

    Returns:
        A tuple of the rendered template, and the token.
    """
    token = get_encoded_token(time_limit=600, arg_dict={"reset_password": user.id})
    return (
        render_template("email/reset_password_text.html", user=user, token=token),
        token,
    )


# send reset email test
def send_notification_test(person: models.Person) -> str:
    """
    Send notification email in testing.

    Args:
        Person: Person to send to.

    Returns:
         Rendered notification text template.
    """
    return render_template("email/notification_text.html", user=person.user)


def verify_encoded_token(token: str, identifier: str) -> models.User:
    """
    Verify token link is valid and return user id.

    Args:
        token: Token to verify.
        identifier: Dictionary key to query after token is decoded. Currently supports "reset_password" or "verify_email"

    Returns:
         User the token identifies.
    """
    try:
        user_id = jwt.decode(
            token, current_app.config["SECRET_KEY"], algorithms=["HS256"]
        )[identifier]
    except Exception:
        return
    return models.User.query.get(user_id)


def get_protocol_type() -> str:
    """
    The desired protocol type depends on the current deployment environment.
    http for local and https if the app is deployed on a remote server

    Returns:
          the active protocol type
    """
    return "http" if current_app.config["MAIL_USE_LOCAL"] == "local" else "https"
