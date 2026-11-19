"""
Workups API (/api/workups) for the workup and purification cards on the reaction page.

Reactions are identified by the names the page uses:
/api/workups/<workgroup>/<workbook>/<reaction code>.

Access:
- Workbook members can read a reaction's steps and the workbook's templates. Anyone
  else gets a 404, so the routes don't reveal which workbooks or reactions exist.
- Only the reaction's creator can change its steps, and only until the reaction is
  locked (403 otherwise).
- Only a template's creator can delete it (403 otherwise).
- Writes only accept JSON bodies (see blueprints/api/utils.py).

Errors are JSON (APIFlask) with a `message` the front end can show.
"""

from typing import Dict, Tuple

from apiflask import abort
from flask import Response, jsonify, request
from flask_login import login_required
from sources import models, services
from sources.extensions import db
from sources.services.inventory import InventoryValidationError
from sources.services.workup import WorkupNotFoundError
from sources.services.workup_steps import WorkupValidationError

from ..utils import json_body
from . import workups_api_bp

LOCKED_MESSAGE = (
    "Only the reaction's creator can change its workup, "
    "and only before the reaction is locked."
)
REACTION_URL = "/<workgroup_name>/<workbook_name>/<reaction_id>"
TEMPLATES_URL = "/<workgroup_name>/<workbook_name>/templates"


@workups_api_bp.errorhandler(WorkupValidationError)
def _handle_validation_error(error: WorkupValidationError) -> Tuple[Response, int]:
    """
    Bad input from the client: discards any half-made changes and returns a 400
    with the service's message.

    Args:
        error: The validation error raised by a workup service.

    Returns:
        The JSON error response and its status code.
    """
    db.session.rollback()
    return jsonify({"message": str(error)}), 400


@workups_api_bp.errorhandler(WorkupNotFoundError)
def _handle_not_found(error: WorkupNotFoundError) -> Tuple[Response, int]:
    """
    An unknown step or template, or one in another reaction or workbook: discards
    any half-made changes and returns a 404.

    Args:
        error: The not-found error raised by a workup service.

    Returns:
        The JSON error response and its status code.
    """
    db.session.rollback()
    return jsonify({"message": "Not found."}), 404


def _member_workbook(workgroup_name: str, workbook_name: str) -> models.WorkBook:
    """
    Gets a workbook the current user is a member of.

    Args:
        workgroup_name: The workgroup name from the URL.
        workbook_name: The workbook name from the URL.

    Returns:
        The workbook.
    """
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    if workbook is None or not services.workbook.is_member(workbook):
        abort(404)
    return workbook


def _member_reaction(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> models.Reaction:
    """
    Gets a reaction in a workbook the current user is a member of.

    Args:
        workgroup_name: The workgroup name from the URL.
        workbook_name: The workbook name from the URL.
        reaction_id: The reaction code from the URL, e.g. "WB1-001".

    Returns:
        The reaction.
    """
    reaction = services.reaction.get_from_names(
        workgroup_name, workbook_name, reaction_id
    )
    if reaction is None:
        abort(404)
    return reaction


def _can_edit(reaction: models.Reaction) -> bool:
    """
    Whether the current user may change the reaction's workup: its creator, while
    the reaction isn't locked (the same rule as the reaction page's autosave).

    Args:
        reaction: The reaction.

    Returns:
        True if the user may change it.
    """
    return (
        services.reaction.check_current_user_is_creator(reaction)
        and reaction.complete != "complete"
    )


def _editable_reaction(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> models.Reaction:
    """
    Gets a reaction whose workup the current user may change.

    Args:
        workgroup_name: The workgroup name from the URL.
        workbook_name: The workbook name from the URL.
        reaction_id: The reaction code from the URL.

    Returns:
        The reaction.
    """
    reaction = _member_reaction(workgroup_name, workbook_name, reaction_id)
    if not _can_edit(reaction):
        abort(403, LOCKED_MESSAGE)
    return reaction


def _steps_response(reaction: models.Reaction, category: models.WorkupCategory) -> Dict:
    """
    Builds the response for a category's steps.

    Args:
        reaction: The reaction.
        category: Workup or purification.

    Returns:
        The steps and whether the user may change them.
    """
    return {
        "category": category.value,
        "canEdit": _can_edit(reaction),
        "steps": [
            step.to_dict() for step in services.workup.list_steps(reaction, category)
        ],
    }


# --- Steps --------------------------------------------------------------------


@workups_api_bp.route("/step-definitions", methods=["GET"])
@login_required
def get_step_definitions() -> Response:
    """Gets the fields each step type records, keyed by step type."""
    return jsonify(services.workup_steps.definitions_as_dict())


@workups_api_bp.route(REACTION_URL, methods=["GET"])
@login_required
def get_steps(workgroup_name: str, workbook_name: str, reaction_id: str) -> Response:
    """Gets a reaction's steps for one category (`?category=WORKUP|PURIFICATION`)."""
    reaction = _member_reaction(workgroup_name, workbook_name, reaction_id)
    category = services.workup.parse_category(request.args.get("category"))
    return jsonify(_steps_response(reaction, category))


@workups_api_bp.route(f"{REACTION_URL}/steps", methods=["POST"])
@login_required
def add_step(workgroup_name: str, workbook_name: str, reaction_id: str) -> Response:
    """Adds an empty step (`{"category": ..., "stepType": ...}`) to the end."""
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    category = services.workup.parse_category(data.get("category"))
    step_type = services.workup.parse_step_type(data.get("stepType"))
    step = services.workup.add_step(reaction, category, step_type)
    return jsonify(step.to_dict()), 201


@workups_api_bp.route(f"{REACTION_URL}/steps/order", methods=["PUT"])
@login_required
def reorder_steps(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> Response:
    """Reorders a category's steps (`{"category": ..., "stepIds": [...]}`)."""
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    category = services.workup.parse_category(data.get("category"))
    services.workup.reorder_steps(reaction, category, data.get("stepIds"))
    return jsonify(_steps_response(reaction, category))


@workups_api_bp.route(f"{REACTION_URL}/steps/<int:step_id>", methods=["PUT"])
@login_required
def update_step(
    workgroup_name: str, workbook_name: str, reaction_id: str, step_id: int
) -> Response:
    """
    Saves a step, replacing its values:
    `{"data": {...}, "notes": ..., "temperature": ..., "materials": [...]}`.
    """
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    step = services.workup.get_step(reaction, step_id)
    step = services.workup.update_step(
        step,
        data.get("data"),
        data.get("notes"),
        data.get("temperature"),
        data.get("materials"),
    )
    return jsonify(step.to_dict())


@workups_api_bp.route(f"{REACTION_URL}/steps/<int:step_id>", methods=["DELETE"])
@login_required
def delete_step(
    workgroup_name: str, workbook_name: str, reaction_id: str, step_id: int
) -> Response:
    """
    Deletes a step. Takes no body: a cross-site DELETE needs a CORS preflight, which
    this app never allows, so another site can't send it.
    """
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    step = services.workup.get_step(reaction, step_id)
    services.workup.delete_step(step)
    return jsonify({"message": "Step deleted."})


@workups_api_bp.route(f"{REACTION_URL}/apply-template", methods=["POST"])
@login_required
def apply_template(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> Response:
    """Appends empty copies of a workbook template's steps (`{"templateId": ...}`)."""
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    try:
        template_id = int(data.get("templateId"))
    except (TypeError, ValueError):
        raise WorkupNotFoundError()
    template = services.workup.get_template(reaction.workbook, template_id)
    services.workup.apply_template(reaction, template)
    return jsonify(_steps_response(reaction, template.category))


@workups_api_bp.route(f"{REACTION_URL}/templates", methods=["POST"])
@login_required
def save_as_template(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> Response:
    """
    Saves a category's step types and order as a workbook template
    (`{"category": ..., "name": ..., "description": ...}`). Any workbook member can
    save another member's workup as a template, since it copies no values.
    """
    reaction = _member_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    category = services.workup.parse_category(data.get("category"))
    template = services.workup.save_as_template(
        reaction, category, data.get("name"), data.get("description")
    )
    return jsonify(template.to_dict()), 201


# --- Material lookups -----------------------------------------------------------


@workups_api_bp.route("/<workgroup_name>/<workbook_name>/solvents", methods=["GET"])
@login_required
def list_solvents(workgroup_name: str, workbook_name: str) -> Response:
    """
    Gets the workbook's solvent list with CHEM21 flags, and the user's CHEM21
    colours, for the coloured solvent picker.
    """
    workbook = _member_workbook(workgroup_name, workbook_name)
    return jsonify(
        {
            "solvents": services.workup_steps.solvent_options(workbook),
            "colours": services.workup_steps.hazard_colours(),
        }
    )


@workups_api_bp.route("/<workgroup_name>/<workbook_name>/compounds", methods=["GET"])
@login_required
def lookup_compound(workgroup_name: str, workbook_name: str) -> Response:
    """
    Looks a compound up by name or CAS number (`?q=`), in the compound database and
    the user's novel compounds in this workgroup. The same lookup as the inventory.
    """
    workbook = _member_workbook(workgroup_name, workbook_name)
    try:
        result = services.inventory.lookup_compound(
            workbook.WorkGroup, request.args.get("q")
        )
    except InventoryValidationError as error:
        abort(400, str(error))
    return jsonify(result)


# --- Templates ----------------------------------------------------------------


@workups_api_bp.route(TEMPLATES_URL, methods=["GET"])
@login_required
def list_templates(workgroup_name: str, workbook_name: str) -> Response:
    """Lists the workbook's workup templates (`?category=` to filter)."""
    workbook = _member_workbook(workgroup_name, workbook_name)
    category_value = request.args.get("category")
    category = (
        services.workup.parse_category(category_value) if category_value else None
    )
    templates = services.workup.list_templates(workbook, category)
    return jsonify(
        [
            {
                **template.to_dict(),
                "canDelete": services.workup.can_delete_template(template),
            }
            for template in templates
        ]
    )


@workups_api_bp.route(f"{TEMPLATES_URL}/<int:template_id>", methods=["DELETE"])
@login_required
def delete_template(
    workgroup_name: str, workbook_name: str, template_id: int
) -> Response:
    """Deletes a template. Only its creator can; reactions that used it keep their steps."""
    workbook = _member_workbook(workgroup_name, workbook_name)
    template = services.workup.get_template(workbook, template_id)
    if not services.workup.can_delete_template(template):
        abort(403, "Only the person who saved this template can delete it.")
    services.workup.delete_template(template)
    return jsonify({"message": "Template deleted."})
