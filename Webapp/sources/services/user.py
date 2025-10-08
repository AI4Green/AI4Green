from datetime import datetime
from typing import List, Optional

from flask_login import current_user
from sources import models, services
from sources.extensions import db
from sqlalchemy import func


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


def add(
    username: str, email: str, fullname: str, password_data: str, person: models.Person
):
    models.User.create(
        username=username,
        email=email,
        fullname=fullname,
        Person=person,
        password_hash=models.User.set_password(password_data),
        privacy_policy_accepted_on=datetime.now(),
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
    return (
        db.session.query(models.User)
        .filter(func.lower(models.User.email) == user_email.lower())
        .first()
    )


def get_retrosynthesis_key(user_email: str) -> Optional[str]:
    """Get the retrosynthesis key associated with a user's email, if they have one.

    Args:
        user_email (str): The user's email.

    Returns:
        Optional[str]: The user's access key to the retrosynthesis service.
    """
    user = (
        db.session.query(models.User)
        .filter(func.lower(models.User.email) == user_email.lower())
        .first()
    )
    if not user:
        return None
    return user.retrosynthesis_access_key.key
