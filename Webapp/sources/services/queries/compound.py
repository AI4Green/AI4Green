from typing import List

from sources import models
from sources.extensions import db


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
