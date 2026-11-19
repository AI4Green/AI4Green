"""
Template instances API (/api/template-instances): a COSHH template filled in for a
reaction or standalone. Reaction workups have their own routes under /api/workups.

Access (see services.template_instances.can_view / can_edit): only the owner.
Others get a 403. Writes accept JSON bodies only, except saving field values, which
needs multipart for file uploads.
"""

import io
from typing import Any, Dict

from apiflask import abort
from flask import Response, json, jsonify, request, send_file
from flask_login import login_required
from sources import db, models, services

from ..utils import json_body
from . import template_instances_api_bp

NOT_OWNER_MESSAGE = "Only the owner can change this form."


def _get_instance(instance_id: int) -> models.TemplateInstance:
    """
    Gets an instance the current user may see.

    Args:
        instance_id: The instance id.

    Returns:
        The instance.
    """
    instance = db.session.get(models.TemplateInstance, instance_id)
    if instance is None:
        abort(404)
    if not services.template_instances.can_view(instance):
        abort(403)
    return instance


def _get_editable_instance(instance_id: int) -> models.TemplateInstance:
    """
    Gets an instance the current user may change.

    Args:
        instance_id: The instance id.

    Returns:
        The instance.
    """
    instance = _get_instance(instance_id)
    if not services.template_instances.can_edit(instance):
        abort(403, NOT_OWNER_MESSAGE)
    return instance


def _reaction_from_names(data: Dict[str, Any]) -> models.Reaction:
    """
    Gets the reaction named in a request, if the user is a member of its workbook.

    Args:
        data: Request data with workgroupName, workbookName and reactionId.

    Returns:
        The reaction.
    """
    reaction = services.reaction.get_from_names(
        str(data.get("workgroupName", "")),
        str(data.get("workbookName", "")),
        str(data.get("reactionId", "")),
    )
    if reaction is None:
        abort(404)
    return reaction


@template_instances_api_bp.route("/", methods=["POST"])
@login_required
def save_new_instance() -> Response:
    """Creates a COSHH form from a template, for a reaction or standalone."""
    data: Dict[str, Any] = json_body()
    template_type = data.get("templateType")
    template_id = data.get("templateId")
    origin = data.get("origin")

    if template_type == "COSHH":
        reaction = None
        if origin == "reactionConstructor":
            reaction = _reaction_from_names(data)
            if not services.reaction.check_current_user_is_creator(reaction):
                abort(404)
        # todo: enforce no duplicates
        new_instance = services.template_instances.add_new_coshh_instance(
            template_id, reaction
        )
        if new_instance is None:
            abort(404)
        return jsonify(new_instance)

    abort(400, "Unknown template type.")


@template_instances_api_bp.route("/<int:instance_id>", methods=["GET"])
@login_required
def get_template_instance(instance_id: int) -> Response:
    """Gets one instance with its sections and field values."""
    instance = _get_instance(instance_id)
    data = instance.to_dict()
    data["canEdit"] = services.template_instances.can_edit(instance)
    return jsonify(data)


@template_instances_api_bp.route("/<int:instance_id>", methods=["PUT"])
@login_required
def save_template_instance(instance_id: int) -> Response:
    """Saves field values (and file uploads) for an instance."""
    instance = _get_editable_instance(instance_id)

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
def download_instance_file(instance_id: int, file_uuid: str) -> Response:
    """Downloads a file saved in one of the instance's fields."""
    instance = _get_instance(instance_id)

    meta = services.field_response.find_file_metadata(instance, file_uuid)
    if meta is None:
        abort(404)

    blob = services.file_attachments.get_blob(file_uuid)
    contents = blob.download_blob().readall()
    if (
        services.file_attachments.sha256_from_file_contents(contents)
        != meta["sha256_checksum"]
    ):
        abort(500, "Integrity check failed")

    return send_file(
        io.BytesIO(contents),
        mimetype=(meta.get("file_details") or {}).get("mimetype"),
        as_attachment=True,
        download_name=meta["display_name"],
    )


@template_instances_api_bp.route("/", methods=["GET"])
@login_required
def get_template_instances() -> Response:
    """Lists instances of one type (?instance_type=COSHH: the current user's forms)."""
    instance_type = request.args.get("instance_type")

    if instance_type == "COSHH":
        return jsonify(services.template_instances.list_coshh_instances_as_dict())

    abort(400, "Unknown instance type.")
