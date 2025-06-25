import json
from typing import Dict, List, Union

from sources import models, services
from sources.extensions import db


def add(
    name: str,
    set_id: str,
    creator: models.Person,
    workbook: models.WorkBook,
    workgroup: models.WorkGroup,
    reactions: List[models.Reaction],
    reactor_dimensions: Dict[str, Union[int, str]],
):
    reaction_set = models.ReactionSet(
        name=name,
        set_id=set_id,
        creator_id=creator.id,
        workbook_id=workbook.id,
        workgroup_id=workgroup.id,
        status="active",
        complete="not complete",
        reactions=reactions,
        reactor_dimensions=json.dumps(reactor_dimensions),
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


def list_active_in_workbook(
    workbook: str, workgroup: str, sort_crit: str = "AZ"
) -> List[models.Reaction]:
    """
    Gets the active reaction_sets for a workbook. Active means the set has not been deleted/archived.

    Args:
        workbook (str): The name of the workbook.
        workgroup (str): The name of the workgroup.
        sort_crit (str, optional): The sorting criteria for the reaction list.
            Defaults to 'AZ' for alphabetical sorting.

    Returns:
        List[models.Reaction]: A list of active reaction_sets in the specified workbook and workgroup,
        sorted based on the specified criteria ('AZ' for alphabetical, 'time' for time of creation).
    """

    query = (
        db.session.query(models.ReactionSet)
        .filter(models.ReactionSet.status == "active")
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
    )
    if sort_crit == "time":
        set_list = query.order_by(models.ReactionSet.time_of_creation.desc()).all()
    elif sort_crit == "AZ":
        set_list = query.order_by(models.ReactionSet.name.asc()).all()
    return set_list


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
        return workbook_abbreviation + "-SET-001"
    most_recent_set_id = newest_set.set_id
    # take the number on the rhs of the set id, remove the 0s, convert to int, add 1, convert to str, add 0s
    new_set_id_number = str(
        int(most_recent_set_id.split("-")[-1].lstrip("0")) + 1
    ).zfill(3)
    new_set_id = workbook_abbreviation + "-SET-" + new_set_id_number
    return new_set_id


def get_from_id(
    set_id: str, workgroup_name: str, workbook_name: str
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
        .filter(models.ReactionSet.set_id == set_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook_name)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )


def to_dict(reaction_set_list: List[models.ReactionSet]) -> List[Dict]:
    """
    Converts a list of reaction sets to a dictionary. Used to render template: '_saved_reactions.html'

    Args:
        reaction_set_list - list of reactions as objects

    Returns:
        A List of dictionaries with the reaction set data required to render the _saved_reactions.html template
    """

    reaction_sets = []
    for idx, reaction_set in enumerate(reaction_set_list):
        description = reaction_set.description
        if len(description) > 250:
            description = description[0:249] + "..."

        if reaction_set.creator.user:
            creator_email = reaction_set.creator.user.email
            creator_username = reaction_set.creator.user.username
        else:
            creator_email = "unknown"
            creator_username = "a deleted profile"

        reaction_set_details = {
            "html_id": idx + 1,
            "name": reaction_set.name,
            "description": description,
            "time_of_creation": str(reaction_set.time_of_creation),
            "time_of_update": str(reaction_set.time_of_update),
            "workgroup": reaction_set.workgroup.name,
            "workbook": reaction_set.workbook.name,
            "completion_status": reaction_set.complete,
            "reaction_id": reaction_set.set_id,
            "creator_email": creator_email,
            "creator_username": creator_username,
            "reactions": services.reaction.to_dict(reaction_set.reactions),
        }
        reaction_sets.append(reaction_set_details)
    return reaction_sets
