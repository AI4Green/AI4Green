from typing import List

from sources import models, services
from sources.extensions import db


def add(
    name: str,
    set_id: str,
    creator: models.Person,
    workbook: models.WorkBook,
    reactions: List[models.Reaction],
):
    reaction_set = models.ReactionSet(
        name=name,
        set_id=set_id,
        creator_id=creator.id,
        workbook_id=workbook.id,
        status="active",
        complete="not complete",
        reactions=reactions,
    )
    db.session.add(reaction_set)
    db.session.commit()
    return reaction_set


def most_recent_in_workbook(workbook_id):
    """
    Retrieves the most recent reaction set in a given workbook.

    Args:
        workbook_id (int): The ID of the workbook.

    Returns:
        models.ReactionSet: The most recent reaction set object in the workbook.
    """
    return (
        db.session.query(models.ReactionSet)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook_id)
        .order_by(models.ReactionSet.set_id.desc())
        .first()
    )


def next_id_in_workbook(workbook_id):
    """
    Generates the next reaction set ID for a given workbook in format WB1-SET-001

    Args:
        workbook_id (int): The ID of the workbook.

    Returns:
        str: The next reaction ID for the workbook.
    """
    workbook = services.workbook.get(workbook_id)
    workbook_abbreviation = workbook.abbreviation
    # find the newest reaction set and then +1 to the id and return
    newest_set = most_recent_in_workbook(workbook_id)
    if not newest_set:
        # if no reactions in workbook yet, then start with 001
        return workbook_abbreviation + "SET-001"
    most_recent_reaction_id = newest_set.reaction_id
    # take the number on the rhs of the set id, remove the 0s, convert to int, add 1, convert to str, add 0s
    new_reaction_id_number = str(
        int(most_recent_reaction_id.split("-")[-1].lstrip("0")) + 1
    ).zfill(3)
    new_set_id = workbook_abbreviation + "-" + new_reaction_id_number
    return new_set_id


def get_from_names(
    set_name: str, workgroup_name: str, workbook_name: str
) -> models.ReactionSet:
    """
    Gets a reaction set with the given name that belongs to the given workgroup and workbook.

    Args:
        set_name (str): The name of the reaction set.
        workgroup_name (str): The name of the workgroup the set belongs to.
        workbook_name (str): The name of the workbook the set belongs to.

    Returns:
        models.ReactionSet: The reaction set with the given name.
    """
    return (
        db.session.query(models.ReactionSet)
        .filter(models.ReactionSet.name == set_name)
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook_name)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )
