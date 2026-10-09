from datetime import datetime
from typing import Dict, List

import pytz
from flask_login import current_user
from sources import models
from sources.extensions import db


def create_new_coshh_template(name: str, description: str) -> Dict:
    """
    Create a new entry in the COSHHTemplate table using name and description variables
    Args:
        name: str, template name
        description: str, template descr

    Returns:
        Models.COSHHTemplate.to_dict(), Dict, dictionary representation of db object
    """
    new_template = models.COSHHTemplate.create(
        name=name,
        description=description,
        template_type=models.template.TemplateType.COSHH,
        time_of_creation=datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        ),
        creator_id=current_user.id,
        institution_id=1,
    )
    db.session.add(new_template)
    db.session.commit()

    return new_template.to_dict()


def list_coshh() -> List[Dict]:
    """
    Queries COSHHTemplate table and returns all templates for current user
    Returns:
        List[models.COSHHTemplate.to_dict()], dictionary representation of all COSHHTemplate objects
    """
    query = db.session.query(models.COSHHTemplate).filter(
        models.Template.creator_id == current_user.id
    )
    return [x.to_dict() for x in query.all()]
