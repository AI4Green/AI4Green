from typing import List

from sqlalchemy import func
from flask_login import current_user
from sources import models
from sources.extensions import db


def list_all() -> List[models.User]:
    """
    Gets a list of all users in the database

    Returns:
         List of all users
    """
    return (
        db.session.query(models.User)
        .order_by(models.User.time_of_creation.desc())
        .all()
    )


def person_from_current_user() -> models.Person:
    """Returns the Person entry in the database Table corresponding to the email of the current user"""
    return (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .first()
    )


def add(
    username: str, email: str, fullname: str, password_data: str, person: models.Person
):
    models.User.create(
        username=username,
        email=email,
        fullname=fullname,
        Person=person,
        password_hash=models.User.set_password(password_data),
    )


def from_id(user_id: int) -> models.User:
    """
    Gets user from User id
    Args:
        user_id: id of user to search for

    Returns:
        models.User with matching id
    """

    return db.session.query(models.User).filter(models.User.id == user_id).first()


def from_email(user_email: str) -> models.User:
    """
        Gets user from User email
        Args:
            user_email: email of user to search for

        Returns:
            models.User with matching email
        """
    return db.session.query(models.User).filter(func.lower(models.User.email) == user_email.lower()).first()