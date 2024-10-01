from sources import services, models
from sources.extensions import db
from typing import List, Optional, Dict, Union
from flask import current_app
from flask_login import current_user
import os
from datetime import datetime
import pytz

def check_reaction_for_controlled_substances(reaction: models.Reaction) -> Union[List[List[str]], None]:
    """
    Checks the provided reaction against AI4Green's list of controlled substances and adds their usage to the database.

    Args:
        reaction (models.Reaction): The reaction to check for controlled substances.

    Returns:
        List[List[str]]: A list of lists containing the InChI for substances that are in AI4Green's controlled substance list.
                        - List[0]: Reactants
                        - List[1]: Reagents
                        - List[2]: Products
                        - List[3]: Solvents
        None: If no controlled substances are used in the reaction.
    """
    substance_smiles = [
        reaction.reactants,
        reaction.reagents,
        reaction.products,
        reaction.solvent
    ]

    # this list of lists maintains the structure of reagents, reactants, products and solvents incase this is useful later
    checks = [check_smiles_in_controlled_substances(smi_list) for smi_list in substance_smiles]

    # if no compounds found return None
    if all(not sublist for sublist in checks):
        return None

    unique_structures = {
        inchi
        for sublist in checks  # flatten the list of lists
        for inchi in sublist
        if not check_duplicate(reaction, inchi)  # filter out duplicates
    }

    if not unique_structures:
        return None

    for substance in unique_structures:
        add(reaction, substance)

    return checks


def check_smiles_in_controlled_substances(smiles_list: List[str]) -> List[str]:
    """
    Checks which of the smiles in the list provided are in AI4Green's controlled substance list.

    Args:
        smiles_list: (List[str]), list of smiles to check

    Returns:
        List of InChI found in controlled substances.
    """
    query_inchis = [services.all_compounds.smiles_to_inchi(smi) for smi in smiles_list]
    return [inchi for inchi in query_inchis if inchi in current_app.config["CONTROLLED_SUBSTANCES"]]


def check_duplicate(reaction: models.Reaction, inchi: str) -> bool:
    """
    Checks db ControlledSubstanceUsage table for entries containing this reaction and InChI and updates the last_edited
    attribute of the entry if duplicate is found

    Args:
        reaction (models.Reaction): The reaction being checked.
        inchi (str): The InChI to check.

    Returns:
        bool, True if duplicate is found, else False

    """
    duplicate = (
        db.session.query(models.ControlledSubstanceUsage)
        .filter(models.ControlledSubstanceUsage.reactions == reaction.id)
        .filter(models.ControlledSubstanceUsage.controlled_substance_inchi == inchi)
        .first()
    )

    if duplicate:
        duplicate.last_edited = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        db.session.commit()
        return True

    return False


def add(reaction: models.Reaction, inchi: str) -> None:
    """
    Adds instance of controlled substance use to the ControlledSubstanceUsage table.

    Args:
        reaction (models.Reaction): The reaction in which the controlled substance is used
        inchi (str): The InChI of the controlled substance

    """
    compound = services.all_compounds.from_inchi(inchi, reaction.workbook)
    location = current_user.most_recent_login_location
    usage = models.ControlledSubstanceUsage(
        creator=reaction.creator,
        workgroups=reaction.workbook.WorkGroup.id,
        workbooks=reaction.workbook.id,
        reactions=reaction.id,
        controlled_substance_name=compound.name,
        controlled_substance_smiles=compound.smiles,
        controlled_substance_cas=compound.cas,
        controlled_substance_inchi=compound.inchi,
        last_location=location["country"]
    )
    db.session.add(usage)
    db.session.commit()

    if location["country"] in current_app.config["EMBARGOED_COUNTRIES"]:
        services.email.send_controlled_substance_alert(inchi, location, reaction)



def list_all() -> List[models.Reaction]:
    """
    Gets a list of all controlled substance usages. For the admin_dashboard

    Returns:
         List of all controlled substance usages
    """
    return (
        (
            db.session.query(models.ControlledSubstanceUsage)
        )
        .order_by(models.ControlledSubstanceUsage.time_of_creation.desc())
        .all()
    )


def controlled_substance_inchi():
    """Returns list of InChI for controlled chemicals that are stores in controlled_substances_inchi.txt"""
    print("150 MISSING INCHI IN THE CONTROLLED SUBSTANCE LIST")
    with open(
            os.path.join(
                os.path.dirname(
                    os.path.dirname(
                        os.path.abspath(__file__)
                    )
                ), "static", "controlled_substances_inchi.txt"
            ), "r"
    ) as f:
        return set(f.read().splitlines())


def uk_arms_embargoes():
    with open(
            os.path.join(
                os.path.dirname(
                    os.path.dirname(
                        os.path.abspath(__file__)
                    )
                ), "static", "UK_arms_embargoed_countries.txt"
            ), "r"
    ) as f:
        return set(f.read().splitlines())

# Add controlled substance list to configs
current_app.config["CONTROLLED_SUBSTANCES"] = controlled_substance_inchi()

# Add embargoed countries list config
current_app.config["EMBARGOED_COUNTRIES"] = uk_arms_embargoes()
