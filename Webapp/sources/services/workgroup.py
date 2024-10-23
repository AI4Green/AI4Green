import re
from typing import List

from flask_login import current_user
from sources import models
from sources.extensions import db
from sqlalchemy import func


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


def get_user_type(workgroup_name: str, user: models.User) -> str:
    """
    Returns the user type of the current user in the specified workgroup
    Args:
        workgroup_name: str, name of workgroup to search
        user: models.User, the user for which to check membership

    Returns:
        user_type: str, user membership type in specified workgroup or None if user is not in workgroup
    """
    user_type = None

    pi = get_workgroup_pi(workgroup_name)
    if user.email in [user.email for user in pi]:
        user_type = "principal_investigator"

    sr = get_workgroup_sr(workgroup_name)
    if user.email in [user.email for user in sr]:
        user_type = "senior_researcher"

    sm = get_workgroup_sm(workgroup_name)
    if user.email in [user.email for user in sm]:
        user_type = "standard_member"

    return user_type


def verify_wg_name(workgroup_name: str, new_name: str) -> str:
    """
    Verifies a name input is unique, is not empty, and contains only alphanumeric characters.
        the new name is stripped of spaces, and lowered to be entirely lowercase. Likewise, all
        names in database are lowered to check for duplicates.
    Args:
            workgroup_name: str, the current workgroup name
            new_name: str, the user input name which requires verification
    Returns:
            name: str, if passes all checks, the name is returned for use
            feedback: str, if check is flagged, feedback is returned with info
    """
    feedback = None
    name_without_spaces = new_name.replace(" ", "")
    duplicates = (
        db.session.query(models.WorkGroup)
        .filter(
            func.lower(func.replace(models.WorkGroup.name, " ", ""))
            == func.lower(name_without_spaces),
            models.WorkGroup.name != workgroup_name,
        )
        .first()
    )
    if duplicates:
        feedback = "Workgroup name already exists"

    elif (
        not new_name.strip()
    ):  # checks if the string is empty or contains only whitespace
        feedback = "Error: Name cannot be empty."

    elif re.search(
        r"[^a-zA-Z0-9 ]", new_name
    ):  # checks for any non-alphanumeric characters in string.
        feedback = "Error: Name contains invalid symbols."

    return feedback if feedback else new_name
