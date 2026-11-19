"""
Sections API (/api/sections): the sections and fields of a COSHH template, edited in
the template builder canvas.

Access: the user must be allowed to manage the section's template
(services.templates.can_manage). Writes only accept JSON bodies.
"""

from typing import Any, Dict

from apiflask import abort
from flask import Response, jsonify
from flask_login import current_user, login_required
from sources import db, models, services
from sources.services.sections import SectionInUseError

from ..utils import json_body
from . import sections_api_bp


def _get_manageable_template(template_id: Any) -> models.Template:
    """
    Gets a template the current user may change.

    Args:
        template_id: The template id from the request.

    Returns:
        The template.
    """
    try:
        template = db.session.get(models.Template, int(template_id))
    except (TypeError, ValueError):
        abort(400, "templateId is required.")
    if template is None or not services.templates.can_view(template):
        abort(404)
    if not services.templates.can_manage(template):
        abort(403)
    return template


def _get_template_section(section_id: int) -> models.Section:
    """
    Gets a template section whose template the current user may change.

    Args:
        section_id: The section id.

    Returns:
        The section.
    """
    section = db.session.get(models.Section, section_id)
    if section is None or section.template_id is None:
        abort(404)
    _get_manageable_template(section.template_id)
    return section


@sections_api_bp.route("/", methods=["GET"])
@login_required
def get_sections() -> Response:
    """Lists the sections of templates the current user created."""
    query = (
        db.session.query(models.Section)
        .join(models.COSHHTemplate)
        .filter(models.COSHHTemplate.creator_id == current_user.id)
    )

    return jsonify([x.to_dict() for x in query.all()])


@sections_api_bp.route("/", methods=["POST"])
@login_required
def save_new_section() -> Response:
    """Adds a section to a template."""
    data: Dict[str, Any] = json_body()
    template = _get_manageable_template(data.get("templateId"))

    section = services.sections.create_section(
        name=data.get("name") or "New Section",
        sort_order=data.get("sortOrder"),
        template=template,
    )
    return jsonify(section.to_dict())


@sections_api_bp.route("/", methods=["PUT"])
@login_required
def update_section() -> Response:
    """Renames and reorders a template's sections, creating any new ones."""
    data: Dict[str, Any] = json_body()
    template = _get_manageable_template(data.get("templateId"))
    sections_data = data.get("sections", [])
    if not isinstance(sections_data, list):
        abort(400, "Send the sections as a list.")

    updated = services.sections.update_template_sections(template, sections_data)
    return jsonify({"sections": [s.to_dict() for s in updated]}), 200


@sections_api_bp.route("/<int:section_id>", methods=["GET"])
@login_required
def get_section_by_id(section_id: int) -> Response:
    """Gets one template section."""
    return jsonify(_get_template_section(section_id).to_dict())


@sections_api_bp.route("/<int:section_id>/fields", methods=["GET"])
@login_required
def get_section_fields(section_id: int) -> Response:
    """Lists a template section's own fields."""
    section = _get_template_section(section_id)
    fields = sorted(section.fields, key=lambda field: field.sort_order or 0)
    return jsonify([x.to_dict() for x in fields])


@sections_api_bp.route("/<int:section_id>/fields", methods=["POST"])
@login_required
def save_new_fields(section_id: int) -> Response:
    """Replaces a template section's fields with the submitted list."""
    section = _get_template_section(section_id)
    fields = json_body(allow_list=True)
    if not isinstance(fields, list):
        abort(400, "Send the fields as a list.")

    try:
        services.sections.save_fields(section, fields)
    except SectionInUseError as error:
        abort(409, str(error))
    return jsonify({"message": "success"}), 201
    # todo: handle multiple input fields, create new select_field_options?
