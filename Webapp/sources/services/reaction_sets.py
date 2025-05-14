from sources import models, services
from sources.extensions import db


def add(name, set_id, creator, workbook_id, reactions):
    reaction_set = models.Reaction(
        name=name,
        reaction_id=set_id,
        creator=creator.id,
        workbooks=workbook_id,
        status="active",
        complete="not complete",
        reactions=reactions,
    )
    db.session.add(reaction_set)
    db.session.commit()
