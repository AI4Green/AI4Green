import uuid
from typing import Dict

from flask_login import current_user
from sources import models, services
from sources.extensions import db


def add_new_coshh_instance(template_id, reaction) -> Dict:
    # ensure new template is a coshh template
    template = services.templates.get_coshh_template_by_id(template_id)
    new_instance = {}
    reaction_id = reaction.id if reaction else None

    if template is not None:
        new_instance = models.COSHHInstance(
            uuid=str(uuid.uuid4()),
            template=template,
            owner_id=current_user.id,
            reaction_id=reaction_id,
            approver_id=current_user.id,  # todo: change this
        )
        db.session.add(new_instance)
        db.session.commit()

    return new_instance.to_dict()


def list_coshh_as_dict():
    """
    Queries COSHHInstance table and returns all coshh instances for the current user.
    Returns:
        List[models.COSHHInstance.to_dict()], dictionary representation of all COSHHInstances objects
    """

    instances = models.COSHHInstance.query.filter_by(owner_id=current_user.id).all()
    return [instance.to_summary_dict() for instance in instances]
