from typing import List

from sources import models
from sources.extensions import db


def get_all_workgroups() -> List[models.WorkGroup]:
    """
    Gets a list of all workgroups in the database

    Returns:
         List of all workgroups
    """
    return db.session.query(models.WorkGroup).all()
