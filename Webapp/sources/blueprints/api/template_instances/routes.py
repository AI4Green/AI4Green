import uuid

from flask import jsonify, request
from flask_login import current_user
from sources import db, models

from . import template_instances_api_bp


@template_instances_api_bp.route("/", methods=["POST"])
def save_new_instance():
    data = request.get_json()

    # check to make sure current user is owner of reaction

    # enforce no duplicates

    print(data)

    new_instance = models.TemplateInstance.create(
        uuid=str(uuid.uuid4()),
        template_type=data.get("templateType", None),
        template_id=data.get("templateId", None),
        owner_id=current_user.id,
        reaction_id=data.get("reactionId", None),
        approver_id=current_user.id,  # todo: change this
    )

    db.session.add(new_instance)
    db.session.commit()

    return jsonify(new_instance.to_dict())
