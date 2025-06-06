from sqlalchemy import func
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


def from_id(person_id: int) -> models.Person:
    return (
        db.session.query(models.Person)
        .filter(models.Person.id == person_id)
        .first()
    )


def from_email(person_email: str) -> models.Person:
    return (
        db.session.query(models.Person)
        .join(models.User)
        .filter(func.lower(models.User.email) == person_email.lower())
        .first()
    )
