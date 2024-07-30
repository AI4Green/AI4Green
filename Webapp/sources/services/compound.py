from typing import List, Optional

from sources import models, services
from sources.extensions import db
from sqlalchemy import func

"""
For all things associated with the Compound table in the database
"""


def count() -> int:
    """
    Gets the number of compounds in the database

    Returns:
        The number of compounds in the database
    """
    return db.session.query(models.Compound).count()


def get_compound_data_error_reports() -> List[models.CompoundDataErrorReport]:
    """
    Gets a list of compound database error reports in the database

    Returns:
         List of all compound database error reports
    """
    # Add time cut off if list grows too large
    return db.session.query(models.CompoundDataErrorReport).all()


def get_smiles(primary_key: int) -> str:
    """
    Gets the smiles from the compound primary key
    """
    return (
        db.session.query(models.Compound.smiles)
        .filter(models.Compound.id == primary_key)
        .first()
    )[0]


def get(primary_key: int) -> models.Compound:
    """Returns the Compound record that matches the primary key ID"""
    return (
        db.session.query(models.Compound)
        .filter(models.Compound.id == primary_key)
        .first()
    )


def get_compound_from_name(name: str) -> models.Compound:
    """
    Retrieves a compound by name.

    Args:
        name: Compound name.

    Returns:
        Compound model.
    """
    return (
        db.session.query(models.Compound)
        .filter(func.lower(models.Compound.name) == name.lower())
        .first()
    )


def get_compound_from_smiles(smiles: str) -> Optional[models.Compound]:
    """
    Retrieve a compound from the database by converting SMILES to InChI.

    Args:
        smiles (str): The SMILES representation of the compound.

    Returns:
        models.Compound: A Compound object corresponding to the provided SMILES.
    """
    inchi = services.all_compounds.smiles_to_inchi(smiles)
    if not inchi:
        return None
    return get_compound_from_inchi(inchi)


def get_compound_from_cas(cas: str) -> models.Compound:
    """
    Retrieves a compound by CAS.

    Args:
        cas: CAS number.

    Returns:
        Compound model.
    """
    return db.session.query(models.Compound).filter(models.Compound.cas == cas).first()


def get_compound_from_inchi(inchi: str) -> models.Compound:
    """
    Retrieve a compound from the database based on its InChI.

    Args:
        inchi (str): The InChI (International Chemical Identifier) of the compound.

    Returns:
        models.Compound: The compound object retrieved from the database.
                        Returns None if no matching compound is found.
    """
    return (
        db.session.query(models.Compound).filter(models.Compound.inchi == inchi).first()
    )
