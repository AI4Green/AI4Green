import uuid
from typing import Dict, List, Optional

from flask_login import current_user
from sources import models, services
from sources.extensions import db


def add_new_coshh_instance(
    template_id: int, reaction: Optional[models.Reaction]
) -> Optional[Dict]:
    """
    Creates a COSHH form from a COSHH template, optionally attached to a reaction.

    Args:
        template_id: The COSHH template id.
        reaction: The reaction to attach it to, or None for a standalone form.

    Returns:
        The new instance as a dict, or None if the template isn't a COSHH template.
    """
    template = services.templates.get_coshh_template_by_id(template_id)
    if template is None:
        return None

    new_instance = models.COSHHInstance(
        uuid=str(uuid.uuid4()),
        template=template,
        owner_id=current_user.id,
        reaction_id=reaction.id if reaction else None,
        approver_id=current_user.id,  # todo: change this
    )
    db.session.add(new_instance)
    db.session.commit()

    return new_instance.to_dict()


def list_coshh_instances_as_dict() -> List[Dict]:
    """
    Queries COSHHInstance table and returns all coshh instances for the current user.
    Returns:
        List[models.COSHHInstance.to_summary_dict()], summaries of the user's COSHH forms
    """

    instances = models.COSHHInstance.query.filter_by(owner_id=current_user.id).all()
    return [instance.to_summary_dict() for instance in instances]


def can_view(instance: models.TemplateInstance) -> bool:
    """
    Whether the current user may see an instance (its owner).

    Args:
        instance: The instance.

    Returns:
        True if the user may see it.
    """
    return instance.owner_id == current_user.id


def can_edit(instance: models.TemplateInstance) -> bool:
    """
    Whether the current user may change an instance (its owner).

    Args:
        instance: The instance.

    Returns:
        True if the user may change it.
    """
    return instance.owner_id == current_user.id
