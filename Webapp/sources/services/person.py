from flask_login import current_user
from sources import models
from sources.extensions import db


def from_current_user_email() -> models.Person:
    """
    Retrieve a Person object associated with the current user's email.

    Returns:
        models.Person: A Person object associated with the current user's email.
    """

    return (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .first()
    )
