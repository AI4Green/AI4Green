from datetime import datetime

import pytz
from flask import jsonify, request
from flask_login import current_user
from sources import db, models

from . import templates_api_bp


@templates_api_bp.route("/", methods=["GET"])
def get_templates():
    template_type = request.args.get("type")
    # todo: move db queries to services?
    query = (
        db.session.query(models.Template)
        .filter(models.Template.creator_id == current_user.id)
        .all()
    )

    if template_type:
        query = query.filter(models.Template.template_type == template_type).all()

    if not query:
        return []

    return jsonify([x.to_dict() for x in query])


@templates_api_bp.route("/<int:template_id>", methods=["GET"])
def get_template(template_id):
    template = models.Template.query.get(template_id)
    return jsonify(template.to_dict())


@templates_api_bp.route("/", methods=["POST"])
def create_new_template():
    data = request.get_json()
    source_id = data.get("source_id", None)
    name = data.get("name", None)
    description = data.get("description", None)

    # todo: handle errors if name or desc are missing
    # todo: include institution id per user

    # if no source id, create a blank template with default values
    if not source_id:
        new_template = models.Template.create(
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

        return jsonify(new_template.to_dict()), 200
