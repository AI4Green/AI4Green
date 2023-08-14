from time import time
from typing import Tuple

import jwt
from flask import current_app, render_template
from sources import models
from sources.extensions import mail


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
    mail.send_email(
        "AI4Chem Reset Your Password",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/reset_password_text.txt", user=user, token=token
        ),
        html_body=render_template(
            "email/reset_password_text.html", user=user, token=token
        ),
    )


def send_notification(person: models.Person):
    """
    Send notifications email to a user.

    Args:
        Person: Person to send to.

    """
    mail.send_email(
        "You have a new AI4Green notification",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[person.user.email],
        text_body=render_template("email/notification_text.txt", user=person.user),
        html_body=render_template("email/notification_text.html", user=person.user),
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
