import re
from datetime import datetime, timedelta
from typing import Dict, List

import pytz
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from sources import models
from sources.extensions import db
from sources.services.ions import (
    reaction_from_ionic_cx_smiles,
    reaction_from_ionic_smiles,
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


def get_from_reaction_id_and_workbook(
    reaction_id: str, workbook_id: int
) -> models.Reaction:
    """
    Gets the reaction from the reaction_id and workbook id

    Args:
        reaction_id - in format WB1-001
        workbook_id - the primary key integer id for the workbook the reaction belongs to
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
    Gets the active reactions for a workbook.

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
            # first we see if it is from marvin js and contains ions
            if len(reaction_smiles.split(" |")) > 1:
                rxn = reaction_from_ionic_cx_smiles(reaction_smiles)
            elif "+" in reaction_smiles or "-" in reaction_smiles:
                rxn = reaction_from_ionic_smiles(reaction_smiles)
                # reactions with no ions - make rxn object directly from string
            else:
                rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
            if size == "small":
                d2d = rdMolDraw2D.MolDraw2DSVG(400, 150)
            else:
                d2d = rdMolDraw2D.MolDraw2DSVG(600, 225)
            d2d.DrawReaction(rxn)
            # return drawing text
            scheme = d2d.GetDrawingText()
            scheme_list.append(scheme)
        else:
            scheme_list.append("")
    return scheme_list


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
