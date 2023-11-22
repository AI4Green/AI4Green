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
