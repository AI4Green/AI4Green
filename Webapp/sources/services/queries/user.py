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
