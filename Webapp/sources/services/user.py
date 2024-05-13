from typing import List

from sources import models
from sources.extensions import db


def list_all() -> List[models.User]:
    """
    Gets a list of all users in the database

    Returns:
         List of all users
    """
    return db.session.query(models.User).all()


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

    return (
        db.session.query(models.User)
        .filter(models.User.id == user_id)
        .first()
    )
