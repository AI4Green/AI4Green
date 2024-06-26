from typing import List

from flask_login import current_user
from sources import models
from sources.extensions import db


def list_all() -> List[models.WorkGroup]:
    """
    Gets a list of all workgroups in the database

    Returns:
         List of all workgroups
    """
    return db.session.query(models.WorkGroup).all()


def from_name(name: str) -> models.WorkGroup:
    """Returns the database workgroup object from the name"""
    return (
        db.session.query(models.WorkGroup).filter(models.WorkGroup.name == name).first()
    )


def get_institution() -> models.Institution:
    """Returns the active institution"""
    # currently all users share institution. May change in future
    return (
        db.session.query(models.Institution)
        .filter(models.Institution.name == "Test User Institution")
        .first()
    )


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


def get_workgroup_from_workgroup_name(workgroup_name: str) -> models.WorkGroup:
    """
    Gets models.WorkGroup object from workgroup name
    Args:
        workgroup_name: str, name of workgroup to return

    Returns:
        models.WorkGroup with matching name
    """
    return (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )


def get_workgroup_pi(workgroup_name: str) -> List[models.User]:
    """
    Gets a list of all Principal Investigator users for the specified workgroup
    Args:
        workgroup_name: name of the workgroup to search

    Returns:
        models.User, list of Users that are Principal Investigators for the specified workgroup
    """

    return (
        db.session.query(models.User.email)
        .join(models.Person)
        .join(models.t_Person_WorkGroup)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .all()
    )


def get_workgroup_sr(workgroup_name: str) -> List[models.User]:
    """
    Gets a list of all Senior Researcher users for the specified workgroup
    Args:
        workgroup_name: name of the workgroup to search

    Returns:
        models.User, list of Users that are Senior Researchers for the specified workgroup
    """

    return (
        db.session.query(models.User)
        .join(models.Person)
        .join(models.t_Person_WorkGroup_2)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .all()
    )


def get_workgroup_sm(workgroup_name: str) -> List[models.User]:
    """
    Gets a list of all Standard Member users for the specified workgroup
    Args:
        workgroup_name: name of the workgroup to search

    Returns:
        models.User, list of Users that are Standard Members for the specified workgroup
    """

    return (
        db.session.query(models.User)
        .join(models.Person)
        .join(models.t_Person_WorkGroup_3)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .all()
    )


def get_user_type(workgroup_name: str) -> str:
    """
    Returns the user type of the current user in the specified workgroup
    Args:
        workgroup_name: str, name of workgroup to search

    Returns:
        user_type: str, user membership type in specified workgroup or None if user is not in workgroup
    """
    user_type = None

    pi = get_workgroup_pi(workgroup_name)
    if current_user.email in [user.email for user in pi]:
        user_type = "principal_investigator"

    sr = get_workgroup_sr(workgroup_name)
    if current_user.email in [user.email for user in sr]:
        user_type = "senior_researcher"

    sm = get_workgroup_sm(workgroup_name)
    if current_user.email in [user.email for user in sm]:
        user_type = "standard_member"

    return user_type
