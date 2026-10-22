import datetime
import json
from typing import Dict, List

from flask import abort
from sources import models
from sources.extensions import db


def serialise(value):
    return value if isinstance(value, str) else json.dumps(value)


def get_or_create_field_response(
    instance: models.TemplateInstance, field: models.Field
) -> models.FieldResponse:
    """
    Finds fields response from provided field and instance, or creates one if it doesnt exist
    Args:
        instance: models.TemplateInstance, this instance to search
        field: models.Field, the field to search

    Returns:
        models.FieldResponse, the matching
    """
    field_response = models.FieldResponse.query.filter_by(
        template_instance_id=instance.id, field_id=field.id
    ).first()
    if field_response is None:
        field_response = models.FieldResponse.create(
            field_id=field.id, template_instance_id=instance.id
        )
        db.session.add(field_response)
        db.session.flush()
    return field_response


def save_value(
    instance: models.TemplateInstance,
    field: models.Field,
    value: [int, str, List, Dict],
) -> models.FieldResponseValue:
    """
    Gets FieldResponse for instance/field combo and saves a FieldResponseValue.

    Args:
        instance:  models.TemplateInstance, the instance to save to
        field: models.Field, the field to associate with the value
        value: int, str, List, Dict, the provided value from front end

    Returns:
        models.FieldResponseValue, the new db entry
    """

    # get or create field response
    field_response = get_or_create_field_response(instance, field)
    serialised = serialise(value)

    # find the latest FieldResponseValue
    latest = (
        models.FieldResponseValue.query.filter_by(field_response_id=field_response.id)
        .order_by(
            models.FieldResponseValue.time_of_response.desc(),
            models.FieldResponseValue.id.desc(),
        )
        .first()
    )
    if latest is not None and latest.value == serialised:
        return None  # unchanged, don't append a duplicate

    # create and add new value to db
    field_response_value = models.FieldResponseValue.create(
        value=serialised,
        field_response_id=field_response.id,
        time_of_response=datetime.now(),
    )
    db.session.add(field_response_value)
    return field_response_value
