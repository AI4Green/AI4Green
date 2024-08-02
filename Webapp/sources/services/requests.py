from typing import List

from flask_login import current_user
from sources import models
from sources.extensions import db


def get_from_email_and_workgroup(email: str, workgroup: models.WorkGroup) -> List[models.WGStatusRequest]:
    """
    Gets requests for user email in provided workgroup
    Args:
        email: str, users email to search
        workgroup: models.WorkGroup, workgroup to search

    Returns:
        models.WGStatusRequest, any requests for that user in the workgroup
    """
    return (
        db.session.query(models.WGStatusRequest)
        .join(models.Person, models.WGStatusRequest.person == models.Person.id)
        .join(models.User)
        .filter(models.User.email == email)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.id == workgroup.id)
        .all()
    )


def find_workgroup_duplicates_for_user(user: models.User, workgroup: models.WorkGroup) -> List[models.WGStatusRequest]:
    """
    Finds duplicate requests by users in a workgroup. To prevent request spamming
    Args:
        user: models.User, user to search for
        workgroup: models.WorkGroup, workgroup to search for

    Returns:
        List of active models.WGStatusRequests for that user in the workgroup
    """
    return (
            db.session.query(models.WGStatusRequest)
            .join(models.Person, models.WGStatusRequest.person == models.Person.id)
            .join(models.User)
            .filter(models.User.email == user.email)
            .join(models.WorkGroup)
            .filter(models.WorkGroup.id == workgroup.id)
            .filter(models.WGStatusRequest.status == 'active')
            .all()
        )


def get_active_in_workgroup_for_pi(user: models.User, workgroup: models.WorkGroup) -> List[models.WGStatusRequest]:
    """
    Finds active requests for principal investigators in a workgroup
    Args:
        user: models.User, user to search for
        workgroup: models.WorkGroup, workgroup to search for

    Returns:
        List of models.WGStatusRequest
    """
    return (
        db.session.query(models.WGStatusRequest)
        .filter(models.WGStatusRequest.status == "active")
        .join(
            models.Person,
            models.WGStatusRequest.principal_investigator == models.Person.id,
        )
        .join(models.User)
        .filter(models.WGStatusRequest.principal_investigator == user.id)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.id == workgroup.id)
        .all()
    )