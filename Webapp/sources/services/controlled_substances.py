import os
from datetime import datetime
from typing import Dict, List, Union

import pytz
from flask import current_app
from sources import models, services
from sources.extensions import db


def check_reaction_for_controlled_substances(
    reaction: models.Reaction,
) -> Union[List[List[str]], None]:
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
        reaction.solvent,
    ]

    # this list of lists maintains the structure of reagents, reactants, products and solvents incase this is useful later
    checks = [
        check_smiles_in_controlled_substances(smi_list) for smi_list in substance_smiles
    ]

    # if no compounds found return None
    if all(not sublist for sublist in checks):
        return None

    location = services.utils.get_location()
    unique_structures = {
        inchi
        for sublist in checks  # flatten the list of lists
        for inchi in sublist
        if not check_duplicate(reaction, inchi, location)  # filter out duplicates
    }

    if not unique_structures:
        return None

    for substance in unique_structures:
        services.controlled_substances.add(reaction, substance, location)

    return checks


def check_smiles_in_controlled_substances(smiles_list: List[str]) -> List[str]:
    """
    Checks which of the smiles in the list provided are in AI4Green's controlled substance list. Doesnt check polymers.

    Args:
        smiles_list: (List[str]), list of smiles to check

    Returns:
        List of InChI found in controlled substances.
    """
    query_inchis = [
        services.all_compounds.smiles_to_inchi(smi) if type(smi) == str else None
        for smi in smiles_list
    ]
    return [
        inchi
        for inchi in query_inchis
        if inchi in current_app.config["CONTROLLED_SUBSTANCES"]
    ]


def check_duplicate(
    reaction: models.Reaction, inchi: str, location: Dict[str, str]
) -> bool:
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
        duplicate.last_location = location["country"]
        db.session.commit()

        detect_usage_in_embargoed_country(reaction, inchi, location)

        return True

    return False


def add(reaction: models.Reaction, inchi: str, location: Dict[str, str]) -> None:
    """
    Adds instance of controlled substance use to the ControlledSubstanceUsage table.

    Args:
        reaction (models.Reaction): The reaction in which the controlled substance is used
        inchi (str): The InChI of the controlled substance
        location (Dict[str: str]): The location of the controlled substance usage

    """
    compound = services.all_compounds.from_inchi(inchi, reaction.workbook)
    usage = models.ControlledSubstanceUsage(
        creator=reaction.creator,
        workgroups=reaction.workbook.WorkGroup.id,
        workbooks=reaction.workbook.id,
        reactions=reaction.id,
        controlled_substance_name=compound.name,
        controlled_substance_smiles=compound.smiles,
        controlled_substance_cas=compound.cas,
        controlled_substance_inchi=compound.inchi,
        last_location=location["country"],
    )
    db.session.add(usage)
    db.session.commit()

    detect_usage_in_embargoed_country(reaction, inchi, location)


def detect_usage_in_embargoed_country(
    reaction: models.Reaction, inchi: str, location: Dict[str, str]
) -> None:
    """
    Checks whether location is on embargoed country list and sends controlled substance usage alert if it is.
    Args:
         reaction (models.Reaction): The reaction in which the controlled substance is used
        inchi (str): The InChI of the controlled substance
        location (Dict[str: str]): The location of the controlled substance usage

    Returns:
        None
    """
    if location["country"] in current_app.config["EMBARGOED_COUNTRIES"]:
        services.email_services.send_controlled_substance_alert(
            inchi, location, reaction
        )


def list_all() -> List[models.ControlledSubstanceUsage]:
    """
    Gets a list of all controlled substance usages. For the admin_dashboard

    Returns:
         List of all controlled substance usages
    """
    return (
        (db.session.query(models.ControlledSubstanceUsage))
        .order_by(models.ControlledSubstanceUsage.time_of_creation.desc())
        .all()
    )


def controlled_substance_inchi():
    """
    Returns list of InChI for controlled chemicals that are stores in controlled_substances_inchi.txt

    Returns:
        List of InChI found in controlled substances.
    """
    with open(
        os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "static",
            "controlled_substances_inchi.txt",
        ),
        "r",
    ) as f:
        return set(f.read().splitlines())


def uk_arms_embargoes():
    """
    Returns list of countries that have a UK arms embargo, listed in UK_arms_embargoed_countries.txt

    Returns:
        List of countries found in UK_arms_embargoed_countries.txt
    """
    with open(
        os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "static",
            "UK_arms_embargoed_countries.txt",
        ),
        "r",
    ) as f:
        return set(f.read().splitlines())
