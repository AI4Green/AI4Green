from typing import List

from flask_login import current_user
from sources import models
from sources.extensions import db


def get_from_email_and_workgroup(email: str, workgroup: models.WorkGroup) -> List[models.Notification]:
    return (
        db.session.query(models.WGStatusRequest)
        .join(models.Person, models.WGStatusRequest.person == models.Person.id)
        .join(models.User)
        .filter(models.User.email == email)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.id == workgroup.id)
        .all()
    )


def find_workgroup_duplicates(email: str, workgroup: models.WorkGroup) -> List[models.WGStatusRequest]:
    return (
        db.session.query(models.WGStatusRequest)
        .join(models.Person, models.WGStatusRequest.person == models.Person.id)
        .join(models.User)
        .filter(models.User.email == email)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.id == workgroup)
        .filter(models.WGStatusRequest.status == 'active')
        .all()
    )


def get_active_in_workgroup(user, workgroup):
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