import io
import uuid
from datetime import datetime

from flask import abort, json, jsonify, request, send_file
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
    origin = data.get("origin", "")
    new_instance = {}

    if origin == "standalone":
        reaction = None

    if origin == "reactionConstructor":
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            workgroup_name, workbook_name
        )

        # check to make sure current user is owner of reaction
        reaction = services.reaction.get_from_reaction_id_and_workbook_id(
            reaction_id, workbook.id
        )

        if not services.reaction.check_current_user_is_creator(reaction):
            abort(404)

    # todo: enforce no duplicates

    if template_type == "COSHH":
        new_instance = services.template_instances.add_new_coshh_instance(
            template_id, reaction
        )

    return jsonify(new_instance)


@template_instances_api_bp.route("/<int:template_id>", methods=["GET"])
@login_required
def get_template_instance(template_id):
    query = models.TemplateInstance.query.get_or_404(template_id)
    if query.owner_id != current_user.id:
        abort(403)
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
        processor.apply_deletions()
        processor.save_pending_files()
        db.session.commit()
        committed = True
    finally:
        if not committed:
            db.session.rollback()
            processor.delete_uploaded_blobs()

    return jsonify({"message": "success!"}), 200


@template_instances_api_bp.route(
    "/<int:instance_id>/files/<file_uuid>", methods=["GET"]
)
@login_required
def download_instance_file(instance_id, file_uuid):
    instance = db.get_or_404(models.TemplateInstance, instance_id)
    if current_user.id != instance.owner_id:
        abort(403)

    meta = services.field_response.find_file_metadata(instance, file_uuid)
    if meta is None:
        abort(404)

    blob = services.file_attachments.get_blob(file_uuid)
    contents = blob.download_blob().readall()
    if (
        services.file_attachments.sha256_from_file_contents(contents)
        != meta["sha256_checksum"]
    ):
        abort(500, description="Integrity check failed")

    return send_file(
        io.BytesIO(contents),
        mimetype=(meta.get("file_details") or {}).get("mimetype"),
        as_attachment=True,
        download_name=meta["display_name"],
    )


@template_instances_api_bp.route("/", methods=["GET"])
@login_required
def get_template_instances():
    instance_type = request.args.get("instance_type", None)
    instance_list = []

    if instance_type == "COSHH":
        instance_list = services.template_instances.list_coshh_as_dict()

    return jsonify(instance_list)
