from typing import List

from flask_login import current_user
from sources import models
from sources.extensions import db


def list_all() -> List[models.WorkBook]:
    """
    Gets a list of all workbooks. Used in admin dashboard.

    Returns:
         List of all workbooks
    """
    return (db.session.query(models.WorkBook)).all()


def get(primary_key: int) -> models.WorkBook:
    """
    Gets the workbook object from the primary key

    Args:
        primary_key: the integer id

    Returns:
        the workbook with the specified primary key
    """
    return (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.id == primary_key)
        .first()
    )


def get_workbooks_from_user_group_combination(workgroup: str) -> List[models.WorkBook]:
    """
    gets list of workbook objects from user and workbook

    Args:
        workgroup: str, workbook to search
        user: flask_login.current_user, current user to retrieve workbooks for

    Returns:
        list of workbook objects for the specified user/workbook combination
    """
    return (
        db.session.query(models.WorkBook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .join(models.t_Person_WorkBook)
        .join(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .all()
    )


def get_workbook_from_group_book_name_combination(
    workgroup_name: str, workbook_name: str
) -> models.WorkBook:
    """
    Gets the workbook object from the workbook name and the name of the workgroup it belongs to.
    Args:
        workgroup_name: the name of the workgroup that the workbook belongs to
        workbook_name: the name of the workbook that we are looking for

    Returns:
        The workbook which matches the name and group combination
    """

    return (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.name == workbook_name)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )


def get_newest_reaction_in_workbooks(
    workbooks: List[models.WorkBook],
) -> models.Reaction:
    """
    Finds the most recent reaction in specified workbooks.
    Args:
        workbooks: list, contains models.Workbook to search from.

    Returns:
        models.Reaction, most recent reaction in workbook
    """

    return (
        db.session.query(models.Reaction)
        .join(models.WorkBook)
        .filter(models.Reaction.workbooks.in_([x.id for x in workbooks]))
        .order_by(models.Reaction.time_of_creation.desc())
        .first()
    )


def get_next_reaction_id_in_workbook(workbook: models.WorkBook) -> str:
    """
    Gets the Id for the next reaction to be added to the workbook
    Args:
        workbook: models.Workbook, workbook reaction will be added to

    Returns:
        Reaction ID: Next reaction ID in the form: workbook_abbreviation-reaction_number
    """
    workbook_abbreviation = workbook.abbreviation
    # find the newest reaction and then +1 to the id and return
    newest_reaction = get_newest_reaction_in_workbooks([workbook])

    if not newest_reaction:
        # if no reactions in workbook yet, then start with 001
        return f"{workbook_abbreviation}-001"
    most_recent_reaction_id = newest_reaction.reaction_id
    # take the number on the rhs of the reaction id, remove the 0s, convert to int, add 1, convert to str, add 0s
    new_reaction_id_number = str(
        int(most_recent_reaction_id.split("-")[-1].lstrip("0")) + 1
    ).zfill(3)
    return f"{workbook_abbreviation}-{new_reaction_id_number}"


def workbooks_from_workgroup(workgroup_name: str) -> List[models.WorkBook]:
    """
    Get a list of workbooks for a given workgroup, for the current user.

    Args:
        workgrou_namep: name of the workgroup

    Returns:
        List of workbooks
    """
    workbooks = (
        db.session.query(models.WorkBook)
        .join(models.WorkGroup)
        .join(models.t_Person_WorkBook)
        .join(models.Person)
        .join(models.User)
        .filter(models.WorkGroup.name == workgroup_name)
        .filter(models.User.email == current_user.email)
        .all()
    )
    return [x for x in workbooks if current_user.Person in x.users]
