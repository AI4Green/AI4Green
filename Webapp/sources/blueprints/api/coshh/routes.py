from datetime import datetime

import pytz
from flask import jsonify, request
from flask_login import current_user
from sources import db, models

from . import coshh_api_bp


@coshh_api_bp.route("/templates", methods=["GET"])
def get_all_templates():
    # todo: sort out this route!
    templates = (
        db.session.query(models.Template)
        .filter(models.Template.creator_id == current_user.id)
        .all()
    )

    if not templates:
        return []

    # need to implement to_dict including permissions and stuff
    return jsonify(templates)

    # return jsonify(
    #     [
    #         {"id": "uu", "name": "ggg", "stage": "Draft", "permissions": "CanPublish"},
    #         {"id": "tt", "name": "ddd", "stage": "Draft", "permissions": "CanPublish"},
    #     ]
    # )
    # return (
    #     db.session.query(models.User)
    #     .filter(func.lower(models.User.email) == user_email.lower())
    #     .first()
    # )


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
