from time import time
from typing import Tuple

import jwt
from flask import current_app, render_template
from sources import models
from sources.extensions import mail


def get_data_export_request_token(
    approver: models.User, data_export_request: models.DataExportRequest
) -> str:
    """
    Get token with expiry time of 1 week.

    Args:
        approver: The user corresponding to a data export request approver
        data_export_request: Data export request to encode.

    Returns:
        A token string.
    """
    expires_in = 604800
    return jwt.encode(
        {
            "data_export_request": data_export_request.id,
            "approver": approver.id,
            "exp": time() + expires_in,
        },
        current_app.config["SECRET_KEY"],
        algorithm="HS256",
    )


def send_data_export_approval_request(
    user: models.User, data_export_request: models.DataExportRequest
) -> None:
    """
    Send password reset email to a given user.

    Args:
        user: User to send to.
        data_export_request: request we are asking PI user to accept or deny
    """
    token = get_reset_password_token(data_export_request)
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Chem Reset Your Password",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/reset_password_text.txt", user=user, token=token, protocol=protocol
        ),
        html_body=render_template(
            "email/reset_password_text.html", user=user, token=token, protocol=protocol
        ),
    )


def get_reset_password_token(user: models.User) -> str:
    """
    Get token with expiry time.

    Args:
        user: User to encode with.

    Returns:
        A token string.
    """
    expires_in = 600
    return jwt.encode(
        {"reset_password": user.id, "exp": time() + expires_in},
        current_app.config["SECRET_KEY"],
        algorithm="HS256",
    )


def send_password_reset(user: models.User) -> None:
    """
    Send password reset email to a given user.

    Args:
        user: User to send to.
    """
    token = get_reset_password_token(user)
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Chem Reset Your Password",
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
    token = get_reset_password_token(user)
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


def verify_reset_password_token(token: str) -> models.User:
    """
    Verify token link is valid and return user id.

    Args:
        token: Token to verify.

    Returns:
         User the token identifies.
    """
    try:
        user_id = jwt.decode(
            token, current_app.config["SECRET_KEY"], algorithms=["HS256"]
        )["reset_password"]
    except Exception:
        return
    return models.User.query.get(user_id)


def verify_data_export_request_token(
    token: str,
) -> Tuple[models.User, models.DataExportRequest]:
    """
    Verify token link is valid and return user id.

    Args:
        token: Token to verify.

    Returns:
         User the token identifies.
    """
    try:
        data_export_request_id, approver_user_id = jwt.decode(
            token, current_app.config["SECRET_KEY"], algorithms=["HS256"]
        )["data_export_request"]["approver"]
        # data_export_request = jwt.decode(token, current_app.config["SECRET_KEY"], algorithms=["HS256"],
        # )
    except Exception:
        return
    return models.User.query.get(approver_user_id), models.DataExportRequest.query.get(
        data_export_request_id
    )


def get_protocol_type() -> str:
    """
    The desired protocol type depends on the current deployment environment.
    http for local and https if the app is deployed on a remote server

    Returns:
          the active protocol type
    """
    return "http" if current_app.config["MAIL_USE_LOCAL"] == "local" else "https"
