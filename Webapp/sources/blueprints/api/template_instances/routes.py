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
    instance = db.get_or_404(models.TemplateInstance, instance_id)

    # only the owner can edit
    if current_user.id != instance.owner_id:
        abort(403)

    data = request.form

    # validate files first
    processor = services.file_attachments.UploadInstanceFiles(request, instance)
    processor.validate_react_submission()

    # Parse JSON strings from form data
    field_responses = json.loads(data.get("fieldResponses", "[]"))
    new_field_responses = json.loads(data.get("newFieldResponses", "[]"))

    # new responses: nfr["id"] is a FIELD id
    for nfr in new_field_responses:
        field = services.field.get_field(instance, nfr.get("id"))
        services.field_response.save_value(instance, field, nfr.get("value"))

    # existing responses: fr["id"] is a FIELD RESPONSE id
    for fr in field_responses:
        field = services.field.get_field_from_field_response(instance, fr.get("id"))
        services.field_response.save_value(instance, field, fr.get("value"))

    # track file save status to catch any errors caused by saving files
    committed = False
    try:
        processor.save_pending_files()
        db.session.commit()
        committed = True
    finally:
        if not committed:
            db.session.rollback()
            processor.delete_uploaded_blobs()

    return jsonify({"message": "success!"}), 200
