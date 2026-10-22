from flask import abort
from sources import models
from sources.extensions import db


def get_field(instance: models.TemplateInstance, field_id: int) -> models.Field:
    """Resolve a Field id, ensuring it belongs to the instance's template."""
    field = db.session.get(models.Field, field_id)
    if field is None or field.section.template_id != instance.template_id:
        abort(400, description="Field does not belong to this template")
    return field


def get_field_from_field_response(
    instance: models.TemplateInstance, response_id: int
) -> models.Field:
    """Resolve a FieldResponse id on this instance to its Field."""
    fr = db.session.get(models.FieldResponse, response_id)
    if fr is None or fr.template_instance_id != instance.id:
        abort(403)
    return get_field(instance, fr.field_id)
