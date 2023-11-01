from typing import List

from sources import models
from sources.extensions import db


def list_all() -> List[models.WorkGroup]:
    """
    Gets a list of all workgroups in the database

    Returns:
         List of all workgroups
    """
    return db.session.query(models.WorkGroup).all()


def get_new_workgroup_requests() -> List[models.WorkGroupRequest]:
    """
    Gets a list of all workgroup requests in the database

    Returns:
         List of all workgroup requests
    """
    return (
        db.session.query(models.WorkGroupRequest)
        .filter(models.WorkGroupRequest.status == "active")
        .all()
    )
