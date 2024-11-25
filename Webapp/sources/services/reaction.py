import re
from datetime import datetime, timedelta
from typing import Dict, List

import pytz
from flask import request
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook
from sources.extensions import db
from sqlalchemy import func


def get_from_name_and_workbook_id(name: str, workbook_id: int) -> models.Reaction:
    """
    Retrieves a reaction based on its name and workbook ID.

    Args:
        name (str): The name of the reaction.
        workbook_id (int): The ID of the workbook to which the reaction belongs.

    Returns:
        models.Reaction: The reaction object matching the name and workbook ID.
    """
    return (
        db.session.query(models.Reaction)
        .filter(func.lower(models.Reaction.name) == name.lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook_id)
        .first()
    )


def get_current_from_request() -> models.Reaction:
    """
    Gets the current reaction for a request from the frontend. Either using request.form or request.json
    Returns:
        The reaction which matches the details of the request
    """
    if request.form:
        reaction = get_current_from_request_form()
    elif request.json:
        reaction = get_current_from_request_json()
    return reaction


def get_current_from_request_form() -> models.Reaction:
    """
    Gets the current reaction using the request.form variable
    Returns:
        Reaction that corresponds to data in request.form
    """
    reaction_id = str(request.form["reactionID"])
    workgroup_name = str(request.form["workgroup"])
    workbook_name = str(request.form["workbook"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    return (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )


def get_current_from_request_json() -> models.Reaction:
    """
    Gets the current reaction using request.json
    Returns:
        Reaction that corresponds to data in request.json
    """
    reaction_id = str(request.json["reactionID"])
    workgroup_name = str(request.json["workgroup"])
    workbook_name = str(request.json["workbook"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    return (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )


def list_recent() -> List[models.Reaction]:
    """
    Gets a list of reactions created in the past 28 days. For the admin_dashboard

    Returns:
         List of all reactions from past 28 days
    """
    cut_off_date = datetime.now(pytz.timezone("Europe/London")).replace(
        tzinfo=None
    ) - timedelta(days=28)
    return (
        (
            db.session.query(models.Reaction).filter(
                models.Reaction.time_of_creation > cut_off_date
            )
        )
        .order_by(models.Reaction.time_of_creation.desc())
        .all()
    )


def count() -> int:
    """
    Gets the number of reactions in the database

    Returns:
        The number of reactions in the database
    """
    return db.session.query(models.Reaction).count()


def get_from_reaction_id_and_workbook_id(
    reaction_id: str, workbook_id: int
) -> models.Reaction:
    """
    Gets the reaction from the reaction_id and workbook id

    Args:
        reaction_id - in format WB1-001
        workbook - The workbook the reaction belongs to

    Returns:
        models.Reaction: The reaction object matching the reaction ID and workbook ID.
    """
    return (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook_id)
        .first()
    )


def list_active_in_workbook(
    workbook: str, workgroup: str, sort_crit: str = "AZ"
) -> List[models.Reaction]:
    """
    Gets the active reactions for a workbook. Active means the reaction has not been deleted/archived.

    Args:
        workbook (str): The name of the workbook.
        workgroup (str): The name of the workgroup.
        sort_crit (str, optional): The sorting criteria for the reaction list.
            Defaults to 'AZ' for alphabetical sorting.

    Returns:
        List[models.Reaction]: A list of active reactions in the specified workbook and workgroup,
        sorted based on the specified criteria ('AZ' for alphabetical, 'time' for time of creation).
    """

    query = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.status == "active")
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
    )
    if sort_crit == "time":
        reaction_list = query.order_by(models.Reaction.time_of_creation.desc()).all()
    elif sort_crit == "AZ":
        reaction_list = query.order_by(models.Reaction.name.asc()).all()
    return reaction_list


def make_scheme_list(reaction_list: List[models.Reaction], size: str) -> List[str]:
    """
    Makes a list of reaction schemes from a list of reactions using the reaction SMILES
    Args:
        reaction_list: list of Reactions objects that we are making scheme images for
        size: the size of scheme image we are making.

    Returns:
        A list of reaction scheme images as strings
    """
    scheme_list = []
    # get rxn schemes
    for reaction in reaction_list:
        reaction_smiles = reaction.reaction_smiles
        # if there are no smiles use a blank string to maintain correct list length
        if re.search("[a-zA-Z]", reaction_smiles):
            # we test to see if ions are present in which case further logic is needed
            scheme_list.append(make_reaction_scheme_image(reaction_smiles, size))
        else:
            scheme_list.append("")
    return scheme_list


def make_reaction_scheme_image(reaction_smiles: str, size: str) -> str:
    """
    Makes a reaction scheme from a reaction smiles according to the size parameter.
    Args:
        reaction_smiles - the SMILES we are making an image of.
        size: whether the image is small or not.
    Returns:
        the image of the reaction scheme as an svg

    """
    # first we see if it is from marvin js and contains ions
    if len(reaction_smiles.split(" |")) > 1:
        rxn = services.ions.reaction_from_ionic_cx_smiles(reaction_smiles)
    elif "+" in reaction_smiles or "-" in reaction_smiles:
        rxn = services.ions.reaction_from_ionic_smiles(reaction_smiles)
        # reactions with no ions - make rxn object directly from string
    else:
        rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
    if size == "small":
        d2d = rdMolDraw2D.MolDraw2DSVG(400, 150)
    else:
        d2d = rdMolDraw2D.MolDraw2DSVG(600, 225)
    d2d.DrawReaction(rxn)
    # return drawing text
    return d2d.GetDrawingText()


def to_dict(reaction_list: List[models.Reaction]) -> List[Dict]:
    """
    Converts a reaction list to a dictionary used to render template: '_saved_reactions.html'

    Args:
        reaction_list - list of reactions as objects
        sort_crit - criteria to sort reactions.

    Returns:
        A List of dictionaries with the reaction data required to render the _saved_reactions.html template

    """

    reactions = []
    for idx, reaction in enumerate(reaction_list):
        # for each reaction get the relevant info and shorten description if it's long
        description = reaction.description
        if reaction.creator_person.user:
            creator_email = reaction.creator_person.user.email
            creator_username = reaction.creator_person.user.username
        else:
            creator_email = "unknown"
            creator_username = "a deleted profile"
        if len(description) > 250:
            description = description[0:249] + "..."
        reaction_details = {
            "html_id": idx + 1,
            "name": reaction.name,
            "description": description,
            "time_of_creation": str(reaction.time_of_creation),
            "time_of_update": str(reaction.time_of_update),
            "reaction_smiles": reaction.reaction_smiles,
            "reaction_table_data": reaction.reaction_table_data,
            "summary_table_data": reaction.summary_table_data,
            "workgroup": reaction.WorkBook.WorkGroup.name,
            "workbook": reaction.WorkBook.name,
            "completion_status": reaction.complete,
            "reaction_id": reaction.reaction_id,
            "creator_email": creator_email,
            "creator_username": creator_username,
            "addenda": reaction.addenda,
        }
        reactions.append(reaction_details)
    return reactions


def add_addendum(
    reaction: models.Reaction, reaction_note_text: str, author: models.Person
) -> models.ReactionNote:
    """
    Add a new addendum to a reaction.

    Args:
        reaction - The reaction object to which the addendum is added.
        reaction_note_text - The text content of the addendum.
        author - The Person object for the author of the addendum

    Returns:
        The newly created reactionNote

    """
    new_addendum = models.ReactionNote(
        text=reaction_note_text,
        time_of_creation=datetime.now(),
        author=author.id,
        reaction=reaction.id,
    )
    db.session.add(new_addendum)
    db.session.commit()
    return new_addendum


def most_recent_in_workbook(workbook_id: int) -> models.Reaction:
    """
    Retrieves the most recent reaction in a given workbook.

    Args:
        workbook_id (int): The ID of the workbook.

    Returns:
        models.Reaction: The most recent reaction object in the workbook.
    """
    return (
        db.session.query(models.Reaction)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook_id)
        .order_by(models.Reaction.reaction_id.desc())
        .first()
    )


def get_next_reaction_id_for_workbook(workbook_id: int) -> str:
    """
    Generates the next reaction ID for a given workbook in format WB1-001

    Args:
        workbook_id (int): The ID of the workbook.

    Returns:
        str: The next reaction ID for the workbook.
    """
    workbook_obj = services.workbook.get(workbook_id)
    workbook_abbreviation = workbook_obj.abbreviation
    # find the newest reaction and then +1 to the id and return
    newest_reaction = most_recent_in_workbook(workbook_id)
    if not newest_reaction:
        # if no reactions in workbook yet, then start with 001
        return workbook_abbreviation + "-001"
    most_recent_reaction_id = newest_reaction.reaction_id
    # take the number on the rhs of the reaction id, remove the 0s, convert to int, add 1, convert to str, add 0s
    new_reaction_id_number = str(
        int(most_recent_reaction_id.split("-")[-1].lstrip("0")) + 1
    ).zfill(3)
    new_reaction_id = workbook_abbreviation + "-" + new_reaction_id_number
    return new_reaction_id


def add(
    name: str,
    reaction_id: str,
    creator: models.Person,
    workbook_id: int,
    reaction_table: Dict[str, any],
    summary_table: Dict[str, any],
    reaction_smiles: str = "",
) -> models.Reaction:
    """
    Adds a reaction to the database.

    Args:
        name (str): The name of the reaction.
        reaction_id (str): The generated id for the reaction in format WB1-001
        creator (models.Person): The creator of the reaction.
        workbook_id (int): The ID of the workbook to which the reaction belongs.
        reaction_table (Dict[str, any]): Data for the reaction table.
        summary_table (Dict[str, any]): Data for the summary table.
        reaction_smiles Optional(str): The SMILES representation of the reaction.

    Returns:
        models.Reaction: The newly added reaction object.
    """
    reaction = models.Reaction(
        name=name,
        reaction_id=reaction_id,
        creator=creator.id,
        workbooks=workbook_id,
        status="active",
        complete="not complete",
        reaction_smiles=reaction_smiles,
        reaction_table_data=reaction_table,
        summary_table_data=summary_table,
    )
    db.session.add(reaction)
    db.session.commit()
    return reaction
