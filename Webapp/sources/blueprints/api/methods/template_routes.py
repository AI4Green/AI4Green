"""
Reaction method templates API, under /api/methods (see
services.reaction_method_templates).

Set this as a separate route, feels maybe more organised than workups?

Access:
- Workbook members can list a workbook's templates and save any reaction's method
  as one. Anyone else gets a 404.
- Applying a template changes the reaction's method, so only its creator can, and
  only until the reaction is locked (403 otherwise).
- Only a template's creator can delete it (403 otherwise).
"""

from apiflask import abort
from flask import Response, jsonify
from flask_login import login_required
from sources import services
from sources.services.reaction_method import MethodNotFoundError

from ..reaction_access import member_reaction, member_workbook
from ..utils import json_body
from . import methods_api_bp
from .routes import REACTION_URL, _editable_reaction, method_response

TEMPLATES_URL = "/<workgroup_name>/<workbook_name>/templates"


@methods_api_bp.route(TEMPLATES_URL, methods=["GET"])
@login_required
def list_templates(workgroup_name: str, workbook_name: str) -> Response:
    """Lists the workbook's method templates, newest first."""
    workbook = member_workbook(workgroup_name, workbook_name)
    return jsonify(
        [
            {
                **template.to_dict(),
                "canDelete": services.reaction_method_templates.can_delete_template(
                    template
                ),
            }
            for template in services.reaction_method_templates.list_templates(workbook)
        ]
    )


@methods_api_bp.route(f"{REACTION_URL}/templates", methods=["POST"])
@login_required
def save_as_template(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> Response:
    """
    Saves the reaction's method as a workbook template
    (`{"name": ..., "description": ...}`).
    """
    reaction = member_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    template = services.reaction_method_templates.save_as_template(
        reaction, data.get("name"), data.get("description")
    )
    return jsonify(template.to_dict()), 201


@methods_api_bp.route(f"{REACTION_URL}/apply-template", methods=["POST"])
@login_required
def apply_template(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> Response:
    """
    Replaces the reaction's method with a workbook template's
    (`{"templateId": ...}`). Returns the method and the template components the
    reaction table doesn't have.
    """
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    try:
        template_id = int(data.get("templateId"))
    except (TypeError, ValueError):
        raise MethodNotFoundError()
    template = services.reaction_method_templates.get_template(
        reaction.workbook, template_id
    )
    result = services.reaction_method_templates.apply_template(reaction, template)
    return jsonify({**method_response(reaction), **result})


@methods_api_bp.route(f"{TEMPLATES_URL}/<int:template_id>", methods=["DELETE"])
@login_required
def delete_template(
    workgroup_name: str, workbook_name: str, template_id: int
) -> Response:
    """Deletes a template. Only its creator can; reactions that used it keep their methods."""
    workbook = member_workbook(workgroup_name, workbook_name)
    template = services.reaction_method_templates.get_template(workbook, template_id)
    if not services.reaction_method_templates.can_delete_template(template):
        abort(403, "Only the person who saved this template can delete it.")
    services.reaction_method_templates.delete_template(template)
    return jsonify({"message": "Template deleted."})
