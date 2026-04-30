import uuid
from datetime import datetime

from flask import json, jsonify, request
from flask_login import current_user
from sources import db, models

from . import template_instances_api_bp


@template_instances_api_bp.route("/", methods=["POST"])
def save_new_instance():
    data = request.get_json()

    # check to make sure current user is owner of reaction

    # todo: enforce no duplicates

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


@template_instances_api_bp.route("/<int:template_id>", methods=["GET"])
def get_template_instance(template_id):
    query = models.TemplateInstance.query.get(template_id)
    data = query.to_dict()  # todo: does this need anymore?

    print(data)
    return jsonify(data)


@template_instances_api_bp.route("/<int:template_id>", methods=["PUT"])
def save_template_instance(template_id):
    data = request.form
    print(data)
    # todo: need to handle files and new update fields, but need to handle template instance reload first
    # Parse JSON strings from form data
    # field_responses = json.loads(data.get("fieldResponses", "[]"))
    new_field_responses = json.loads(data.get("newFieldResponses", "[]"))
    # file_field_responses = json.loads(data.get("fileFieldResponses", "[]"))
    # new_file_field_responses = json.loads(data.get("newFileFieldResponses", "[]"))

    # section_id = data.get("sectionId")
    template_instance_id = data.get("recordId")

    for r in new_field_responses:
        print(r)
        field_id = r.get("id")
        value = r.get("value")

        field = models.Field.query.get(field_id)

        field_response = models.FieldResponse.create(
            field_id=field.id,
            template_instance_id=template_instance_id,
        )

        field_response_value = models.FieldResponseValue.create(
            value=value if isinstance(value, str) else json.dumps(value),
            field_response_id=field_response.id,
            time_of_response=datetime.now(),
        )

        db.session.add(field_response, field_response_value)
        db.session.commit()

        return jsonify({"message": "success!"}), 200

    # query = models.TemplateInstance.query.get(template_id)
    # data = query.to_dict() # todo: does this need anymore?
    # return jsonify(data)
