from datetime import datetime
from typing import Dict, List, Optional

import pytz
from flask_login import current_user
from sources import models, services
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


def list_coshh_templates_as_dict() -> List[Dict]:
    """
    Queries COSHHTemplate table and returns all templates. Templates can only be made by H&S admin.
    Returns:
        List[models.COSHHTemplate.to_dict()], dictionary representation of all COSHHTemplate objects
    """

    query = db.session.query(models.COSHHTemplate)
    return [x.to_dict() for x in query.all()]


def get_coshh_template_by_id(template_id: int) -> Optional[models.COSHHTemplate]:
    """
    Get COSHHTemplate by template id
    Args:
        template_id: str, the template id

    Returns:
        models.COSHHTemplate, or None if there is no COSHH template with that id
    """
    return db.session.get(models.COSHHTemplate, template_id)


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


def can_view(template: models.Template) -> bool:
    """
    Whether the current user may see a template.

    - Workup templates: members of the template's workbook.
    - Other (COSHH) templates: H&S admins.

    Args:
        template: The template.

    Returns:
        True if the user may see it.
    """
    if isinstance(template, models.WorkupTemplate):
        return services.workbook.is_member(template.workbook)
    return is_hs_admin()


def can_manage(template: models.Template) -> bool:
    """
    Whether the current user may change a template's sections and fields.

    Same rule as can_view for now; kept separate so it can tighten later
    (e.g. once publishing/staging is built).

    Args:
        template: The template.

    Returns:
        True if the user may change it.
    """
    return can_view(template)
