from typing import List

from sources import models
from sources.extensions import db


def list_all() -> List[models.WorkBook]:
    """
    Gets a list of all workbooks

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


def get_next_reaction_id_in_workbook(workbook: models.WorkBook) -> str:
    workbook_abbreviation = workbook.abbreviation
    # find the newest reaction and then +1 to the id and return
    newest_reaction = (
        db.session.query(models.Reaction)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .order_by(models.Reaction.reaction_id.desc())
        .first()
    )
    if not newest_reaction:
        # if no reactions in workbook yet, then start with 001
        return f"{workbook_abbreviation}-001"
    most_recent_reaction_id = newest_reaction.reaction_id
    # take the number on the rhs of the reaction id, remove the 0s, convert to int, add 1, convert to str, add 0s
    new_reaction_id_number = str(
        int(most_recent_reaction_id.split("-")[-1].lstrip("0")) + 1
    ).zfill(3)
    return f"{workbook_abbreviation}-{new_reaction_id_number}"