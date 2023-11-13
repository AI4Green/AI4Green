from typing import Tuple

from flask import abort, request
from flask_login import current_user
from sources import models
from sources.auxiliary import security_member_workgroup_workbook
from sources.extensions import db
from sources.services import queries


def get_smiles(primary_key: Tuple[str, int]) -> str:
    """
    Gets the novel compound's SMILES string from the primary key if the entry has the SMILES attribute

    Returns:
        The SMILES string corresponding to the primary key or None
    """
    primary_key = (primary_key[0], int(primary_key[1]))
    workbook = queries.workbook.get_workbook_from_primary_key(primary_key[1])
    if current_user.Person not in workbook.users:
        abort(401)

    return (
        db.session.query(models.NovelCompound.smiles)
        .filter(models.NovelCompound.name == primary_key[0])
        .join(models.WorkBook)
        .filter(models.WorkBook.id == primary_key[1])
        .first()
    )[0]
