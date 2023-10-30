from typing import Tuple

from flask import abort
from flask_login import current_user
from sources import models
from sources.extensions import db
from sources.services import queries


def get_smiles_from_primary_key(primary_key: Tuple[str, int]) -> int:
    """
    Gets the novel compound's SMILES string from the primary key if the entry has the SMILES attribute

    Returns:
        The SMILES string corresponding to the primary key or None
    """
    workbook = queries.workbook.get_workbook_from_primary_key(primary_key[0])
    if current_user.email not in workbook.users:
        abort(401)

    return (
        db.session.query(models.NovelCompound.smiles)
        .filter(models.NovelCompound.id == primary_key[0])
        .join(models.WorkBook)
        .filter(models.WorkBook.id == primary_key[1])
        .first()
    )
