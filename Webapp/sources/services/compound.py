from typing import List

from rdkit import Chem
from sources import models
from sources.extensions import db
from sqlalchemy import func


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


def get_compound_from_inchi(inchi: str) -> models.Compound:
    """
    Retrieves a compound by InChI

    Args:
        inchi: compound's InChI (structurally derived identifier)

    Returns:
        compound model.
    """
    return (
        db.session.query(models.Compound).filter(models.Compound.inchi == inchi).first()
    )


def get_compound_from_smiles(smiles: str) -> models.Compound:
    """
    Retrieves a compound from SMILES by converting to InChI via an RDKit mol object

    :param smiles:
    :return:
    """
    mol = Chem.MolFromSmiles(smiles)
    inchi = Chem.MolToInchi(mol)
    return (
        db.session.query(models.Compound).filter(models.Compound.inchi == inchi).first()
    )


def get_compound_from_cas(cas: str) -> models.Compound:
    """
    Retrieves a compound by CAS.

    Args:
        cas: CAS number.

    Returns:
        Compound model.
    """
    return db.session.query(models.Compound).filter(models.Compound.cas == cas).first()
