"""
Templates API (/api/templates) for COSHH templates. Workup templates have their own
routes under /api/workups.

Access:
- COSHH templates are managed by H&S admins (403 for anyone else). Any logged-in
  user can list them, to pick one for a COSHH form.
- Writes only accept JSON bodies (see blueprints/api/utils.py).
"""

from typing import Any, Dict

from apiflask import abort
from flask import Response, jsonify, request
from flask_login import login_required
from sources import models, services
from sources.extensions import db

from ..utils import json_body
from . import templates_api_bp


def _get_template(template_id: int) -> models.Template:
    """
    Gets a template the current user may see.

    Args:
        template_id: The template id.

    Returns:
        The template.
    """
    template = db.session.get(models.Template, template_id)
    if template is None:
        abort(404)
    if not services.templates.can_view(template):
        abort(403)
    return template


@templates_api_bp.route("/", methods=["GET"])
@login_required
def get_templates() -> Response:
    """Lists templates of one type (?template_type=COSHH)."""
    if request.args.get("template_type") == "COSHH":
        return jsonify(services.templates.list_coshh_templates_as_dict())
    abort(400, "Unknown template type.")


@templates_api_bp.route("/<int:template_id>", methods=["GET"])
@login_required
def get_template(template_id: int) -> Response:
    """Gets one template."""
    return jsonify(_get_template(template_id).to_dict())


@templates_api_bp.route("/", methods=["POST"])
@login_required
def create_new_template() -> Response:
    """Creates an empty template of the given type."""
    data: Dict[str, Any] = json_body()
    name = (data.get("name") or "").strip()
    description = data.get("description") or ""
    template_type = data.get("templateType")

    if not name:
        abort(400, "A template needs a name.")

    if template_type == "COSHH":
        if not services.templates.is_hs_admin():
            abort(403, "Only H&S admins can create COSHH templates.")
        return jsonify(services.templates.create_new_coshh_template(name, description))

    abort(400, "Unknown template type.")


@templates_api_bp.route("/<int:template_id>", methods=["PUT"])
@login_required
def edit_template(template_id: int) -> Response:
    """Renames a template or changes its description."""
    template = _get_template(template_id)
    if not services.templates.can_manage(template):
        abort(403)
    data: Dict[str, Any] = json_body()

    name = (data.get("name") or "").strip()
    if not name:
        abort(400, "A template needs a name.")
    template.name = name
    template.description = data.get("description") or ""
    db.session.commit()
    return jsonify(template.to_dict())


@templates_api_bp.route("/<int:template_id>", methods=["DELETE"])
@login_required
def delete_template(template_id: int) -> Response:
    """Deletes a template that hasn't been used yet."""
    template = _get_template(template_id)
    if not services.templates.can_manage(template):
        abort(403)
    if template.template_instances:
        abort(409, "This template is in use, so it can't be deleted.")
    db.session.delete(template)
    db.session.commit()
    return jsonify({"message": "Template deleted."})


@templates_api_bp.route("/<int:template_id>/sections", methods=["GET"])
@login_required
def get_template_sections(template_id: int) -> Response:
    """Lists a template's sections in order."""
    template = _get_template(template_id)
    sections = sorted(template.sections, key=lambda x: (x.sort_order or 0, x.id))
    return jsonify([x.to_dict() for x in sections])
