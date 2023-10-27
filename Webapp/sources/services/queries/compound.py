from sources import models
from sources.extensions import db


def get_compound_count() -> int:
    """
    Gets the number of compounds in the database

    Returns:
        The number of compounds in the database
    """
    return db.session.query(models.Compound).count()


def get_smiles_from_pk(int) -> str:
    """
    Gets the smiles from the compound primary key

    """
