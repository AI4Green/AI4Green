from sources import models, services
from sources.extensions import db


def get_coshh_form_from_reaction_id(reaction_id):
    return (
        db.session.query(models.TemplateInstance)
        .filter(models.TemplateInstance.template_type == "COSHH")
        .join(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .first()
    )
