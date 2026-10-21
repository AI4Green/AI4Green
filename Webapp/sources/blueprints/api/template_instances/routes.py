import uuid
from datetime import datetime

from flask import abort, json, jsonify, request
from flask_login import current_user, login_required
from sources import db, models, services

from . import template_instances_api_bp


@template_instances_api_bp.route("/", methods=["POST"])
@login_required
def save_new_instance():
    data = request.get_json()

    reaction_id = data.get("reactionId", "")
    workgroup_name = data.get("workgroupName", "")
    workbook_name = data.get("workbookName", "")
    template_type = data.get("templateType", "")
    template_id = data.get("templateId", "")
    new_instance = {}

    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )

    # check to make sure current user is owner of reaction
    reaction = services.reaction.get_from_reaction_id_and_workbook_id(
        reaction_id, workbook.id
    )

    # todo: enforce no duplicates

    if template_type == "COSHH":
        new_instance = services.template_instances.add_new_coshh_instance(
            template_id, reaction
        )

    return jsonify(new_instance)


@template_instances_api_bp.route("/<int:template_id>", methods=["GET"])
@login_required
def get_template_instance(template_id):
    query = models.TemplateInstance.query.get(template_id)
    data = query.to_dict()

    return jsonify(data)


@template_instances_api_bp.route("/<int:instance_id>", methods=["PUT"])
@login_required
def save_template_instance(instance_id):
    instance = models.TemplateInstance.query.get_or_404(instance_id)

    # TODO: enforce who may edit this instance, e.g.
    if current_user.id not in (
        instance.owner_id,
        getattr(instance, "approver_id", None),
    ):
        abort(403)

    data = request.form

    # Parse JSON strings from form data
    field_responses = json.loads(data.get("fieldResponses", "[]"))
    new_field_responses = json.loads(data.get("newFieldResponses", "[]"))

    # todo: implement files
    # file_field_responses = json.loads(data.get("fileFieldResponses", "[]"))
    # new_file_field_responses = json.loads(data.get("newFileFieldResponses", "[]"))

    # processor = services.file_attachments.UploadExperimentDataFiles(request)
    # processor.process_react_submission()

    def serialise(value):
        return value if isinstance(value, str) else json.dumps(value)

    def save_value(field_id, value):
        """Append a value to this instance's response for the field,
        creating the response if it doesn't exist yet."""
        field = models.Field.query.get(field_id)
        if field is None or field.section.template_id != instance.template_id:
            abort(400, description="Field does not belong to this template")

        field_response_object = models.FieldResponse.query.filter_by(
            template_instance_id=instance.id, field_id=field.id
        ).first()

        if field_response_object is None:
            field_response_object = models.FieldResponse.create(
                field_id=field.id,
                template_instance_id=instance.id,
            )
            db.session.add(field_response_object)
            db.session.flush()

        serialised = serialise(value)

        latest = (
            models.FieldResponseValue.query.filter_by(
                field_response_id=field_response_object.id
            )
            .order_by(
                models.FieldResponseValue.time_of_response.desc(),
                models.FieldResponseValue.id.desc(),
            )
            .first()
        )
        if latest is not None and latest.value == serialised:
            return  # unchanged, don't append a duplicate

        field_response_value = models.FieldResponseValue.create(
            value=serialised,
            field_response_id=field_response_object.id,
            time_of_response=datetime.now(),
        )
        db.session.add(field_response_value)

    # new responses: r["id"] is a FIELD id
    for r in new_field_responses:
        save_value(r.get("id"), r.get("value"))

    # existing responses: r["id"] is a FIELD RESPONSE id
    for r in field_responses:
        field_response = models.FieldResponse.query.get_or_404(r.get("id"))
        if field_response.template_instance_id != instance.id:
            abort(403)
        save_value(field_response.field_id, r.get("value"))

    db.session.commit()

    return jsonify({"message": "success!"}), 200
