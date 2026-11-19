from datetime import datetime
from typing import Dict, List

import pytz
from flask_login import current_user
from sources import models
from sources.extensions import db

HS_ADMIN_ROLE = "H&S Admin"


def create_new_coshh_template(name: str, description: str) -> Dict:
    """
    Create a new entry in the COSHHTemplate table using name and description variables
    Args:
        name: str, template name
        description: str, template descr

    Returns:
        Models.COSHHTemplate.to_dict(), Dict, dictionary representation of db object
    """
    new_template = models.COSHHTemplate(
        name=name,
        description=description,
        time_of_creation=datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        ),
        creator_id=current_user.id,
        institution_id=1,
    )
    db.session.add(new_template)
    db.session.commit()

    return new_template.to_dict()


def list_coshh_as_dict() -> List[Dict]:
    """
    Queries COSHHTemplate table and returns all templates. Templates can only be made by H&S admin.
    Returns:
        List[models.COSHHTemplate.to_dict()], dictionary representation of all COSHHTemplate objects
    """

    query = db.session.query(models.COSHHTemplate)
    return [x.to_dict() for x in query.all()]


def get_coshh_template_by_id(template_id: int) -> Dict:
    """
    Get COSHHTemplate by template id
    Args:
        template_id: str, the template id

    Returns:
        models.COSHHTemplate
    """
    template = db.session.get(models.COSHHTemplate, template_id)
    return template if template else None


def is_hs_admin() -> bool:
    """
    Whether the current user is an H&S admin, who manages COSHH templates.

    Returns:
        True for H&S admins.
    """
    return (
        current_user.is_authenticated
        and current_user.Role is not None
        and current_user.Role.name == HS_ADMIN_ROLE
    )
