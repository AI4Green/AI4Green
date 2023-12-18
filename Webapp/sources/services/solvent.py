from typing import List

from sources import models
from sources.extensions import db


def add(name: str, hazards: str, nc: models.NovelCompound) -> models.Solvent:
    """
    Creates a solvent model in the database.

    Returns:
        Solvent model.
    """
    solvent = models.Solvent(
        name=name,
        flag=5,
        hazard=hazards,
        novel_compound=[nc],
    )
    db.session.add(solvent)
    db.session.commit()
    return solvent


def get_default_list() -> List[models.Solvent]:
    """
    Gets the default list of solvents. I.e., excluding any in a workbook. Used in the demo mode.

    Returns:
        The default list of solvents from CHEM21
    """
    sol_rows = db.session.query(models.Solvent).all()
    return [x for x in sol_rows if x.novel_compound == []]


def get_workbook_list(workbook: models.WorkBook) -> List[models.Solvent]:
    """
    Gets the list of solvents for the active workbook

    Args:
        workbook object of the active workbook

    Returns:
        The default list of solvents from CHEM21 and any solvents added to the workbook
    """
    sol_rows = db.session.query(models.Solvent).all()
    return [
        x
        for x in sol_rows
        if x.novel_compound == []
        or (workbook is not None and x.novel_compound[0].workbook == workbook.id)
    ]
