from datetime import datetime

import pytz
from flask import jsonify, request
from flask_login import current_user
from sources import db, models

from . import coshh_api_bp


@coshh_api_bp.route("/templates", methods=["GET"])
def get_all_templates():
    templates = (
        db.session.query(models.Template)
        .filter(models.Template.creator_id == current_user.id)
        .all()
    )

    if not templates:
        return []

    return jsonify([x.to_dict() for x in templates])


@coshh_api_bp.route("/templates", methods=["POST"])
def create_new_template():
    data = request.get_json()
    source_id = data.get("source_id", None)
    name = data.get("name", None)
    description = data.get("description", None)

    # todo: handle errors if name or desc are missing
    # todo: include institution id per user

    # if no source id, create a blank template with default values
    if not source_id:
        template = models.Template.create(
            name=name,
            description=description,
            template_type=models.template.TemplateType.COSHH,
            time_of_creation=datetime.now(pytz.timezone("Europe/London")).replace(
                tzinfo=None
            ),
            creator_id=current_user.id,
            institution_id=1,
        )
        db.session.add(template)
        db.session.commit()

        return jsonify({"message": "blank template created"}), 200


@coshh_api_bp.route("/templates/<int:template_id>", methods=["GET"])
def get_template(template_id):
    template = models.Template.query.get(template_id)
    return jsonify(template.to_dict())


@coshh_api_bp.route("/templates/section_types", methods=["GET"])
def get_template_section_types():
    section_types = models.SectionType.query.all()
    return jsonify([x.to_dict for x in section_types])
