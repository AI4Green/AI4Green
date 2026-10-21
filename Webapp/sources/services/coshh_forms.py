from sources import models, services
from sources.extensions import db


def get_coshh_form_from_reaction_pk(reaction_pk):
    return (
        db.session.query(models.COSHHInstance)
        .filter(models.COSHHInstance.reaction_id == reaction_pk)
        .first()
    )
