"""
Reaction methods API (/api/methods) for the method card on the reaction page: the
reaction's setup, conditions and order of addition (see services.reaction_method).

Reactions are identified by the names the page uses:
/api/methods/<workgroup>/<workbook>/<reaction code>.

Access (see blueprints/api/reaction_access.py):
- Workbook members can read a reaction's method. Anyone else gets a 404.
- Only the reaction's creator can change it, and only until the reaction is locked
  (403 otherwise).
- Writes only accept JSON bodies (see blueprints/api/utils.py).

Every write returns the whole method, as the GET does, so the card can show the
saved values. Errors are JSON (APIFlask) with a `message` the front end can show.
"""

from typing import Any, Dict, Tuple

from flask import Response, jsonify
from flask_login import login_required
from sources import models, services
from sources.extensions import db
from sources.services.reaction_method import MethodNotFoundError, MethodValidationError

from ..reaction_access import can_edit, editable_reaction, member_reaction
from ..utils import json_body
from . import methods_api_bp

LOCKED_MESSAGE = (
    "Only the reaction's creator can change its method, "
    "and only before the reaction is locked."
)
REACTION_URL = "/<workgroup_name>/<workbook_name>/<reaction_id>"


# error handlers for workups/methods could probably be condensed, but feels like its clearer like this
@methods_api_bp.errorhandler(MethodValidationError)
def _handle_validation_error(error: MethodValidationError) -> Tuple[Response, int]:
    """
    Bad input from the client: discards any half-made changes and returns a 400
    with the service's message.

    Args:
        error: The validation error raised by a method service.

    Returns:
        The JSON error response and its status code.
    """
    db.session.rollback()
    return jsonify({"message": str(error)}), 400


@methods_api_bp.errorhandler(MethodNotFoundError)
def _handle_not_found(error: MethodNotFoundError) -> Tuple[Response, int]:
    """
    An unknown addition or template, or one in another reaction or workbook:
    discards any half-made changes and returns a 404.

    Args:
        error: The not-found error raised by a method service.

    Returns:
        The JSON error response and its status code.
    """
    db.session.rollback()
    return jsonify({"message": "Not found."}), 404


def _editable_reaction(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> models.Reaction:
    """
    Gets a reaction whose method the current user may change.

    Args:
        workgroup_name: The workgroup name from the URL.
        workbook_name: The workbook name from the URL.
        reaction_id: The reaction code from the URL.

    Returns:
        The reaction.
    """
    return editable_reaction(workgroup_name, workbook_name, reaction_id, LOCKED_MESSAGE)


def method_response(reaction: models.Reaction) -> Dict[str, Any]:
    """
    Builds the response for a reaction's method.

    Args:
        reaction: The reaction.

    Returns:
        The method, its additions, the reaction table's components, and whether
        the user may change them.
    """
    return {
        "canEdit": can_edit(reaction),
        **services.reaction_method.serialise(reaction),
    }


@methods_api_bp.route("/options", methods=["GET"])
@login_required
def get_options() -> Response:
    """Gets the choices for each of the method's dropdowns, with their labels."""
    return jsonify(services.reaction_method.options())


@methods_api_bp.route(REACTION_URL, methods=["GET"])
@login_required
def get_method(workgroup_name: str, workbook_name: str, reaction_id: str) -> Response:
    """Gets a reaction's method, with amounts from its reaction table."""
    reaction = member_reaction(workgroup_name, workbook_name, reaction_id)
    return jsonify(method_response(reaction))


@methods_api_bp.route(f"{REACTION_URL}/procedure-text", methods=["GET"])
@login_required
def get_procedure_text(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> Response:
    """Gets the method written out as experimental procedure text."""
    reaction = member_reaction(workgroup_name, workbook_name, reaction_id)
    method = services.reaction_method.serialise(reaction)
    return jsonify({"text": services.procedure_text.method_text(method)})


@methods_api_bp.route(f"{REACTION_URL}/setup", methods=["PUT"])
@login_required
def update_setup(workgroup_name: str, workbook_name: str, reaction_id: str) -> Response:
    """
    Saves the setup and conditions (`{"setup": {...}, "conditions": {...}}`). Only
    the keys sent are changed.
    """
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    services.reaction_method.update_setup_conditions(
        reaction, data.get("setup"), data.get("conditions")
    )
    return jsonify(method_response(reaction))


@methods_api_bp.route(f"{REACTION_URL}/additions", methods=["POST"])
@login_required
def add_addition(workgroup_name: str, workbook_name: str, reaction_id: str) -> Response:
    """
    Adds an addition of reaction table components to the end
    (`{"components": [{"role": ..., "key": ...}]}`).
    """
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    services.reaction_method.add_addition(reaction, data.get("components"))
    return jsonify(method_response(reaction)), 201


@methods_api_bp.route(f"{REACTION_URL}/additions/order", methods=["PUT"])
@login_required
def reorder_additions(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> Response:
    """Reorders the additions (`{"additionIds": [...]}`)."""
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    services.reaction_method.reorder_additions(reaction, data.get("additionIds"))
    return jsonify(method_response(reaction))


@methods_api_bp.route(f"{REACTION_URL}/additions/<int:addition_id>", methods=["PUT"])
@login_required
def update_addition(
    workgroup_name: str, workbook_name: str, reaction_id: str, addition_id: int
) -> Response:
    """Saves an addition. Only the keys sent are changed."""
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    data = json_body()
    addition = services.reaction_method.get_addition(reaction, addition_id)
    services.reaction_method.update_addition(addition, data)
    return jsonify(method_response(reaction))


@methods_api_bp.route(f"{REACTION_URL}/additions/<int:addition_id>", methods=["DELETE"])
@login_required
def delete_addition(
    workgroup_name: str, workbook_name: str, reaction_id: str, addition_id: int
) -> Response:
    """
    Deletes an addition. Takes no body: a cross-site DELETE needs a CORS preflight,
    which this app never allows, so another site can't send it.
    """
    reaction = _editable_reaction(workgroup_name, workbook_name, reaction_id)
    addition = services.reaction_method.get_addition(reaction, addition_id)
    services.reaction_method.delete_addition(addition)
    return jsonify(method_response(reaction))
