from typing import Tuple

from sources import models
from sources.extensions import db


def get_smiles_from_pk(primary_key: Tuple[str, int]) -> int:
    """
    Gets the novel compound's SMILES string from the primary key if the entry has the SMILES attribute

    Returns:
        The SMILES string corresponding to the primary key
    """

    return db.session.query(models.Compound).count()
