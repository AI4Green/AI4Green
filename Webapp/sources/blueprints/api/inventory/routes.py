"""
Inventory API for the React inventory pages.

Access:
- Every route needs a logged-in member of the workgroup. Anyone else gets a
  404, so the routes don't reveal which workgroups exist.
- Reading is open to all members; changes need a principal investigator or
  senior researcher (services.inventory.EDITOR_ROLES), otherwise 403.
- Writes only accept JSON bodies. App-wide CSRF tokens are disabled, and a
  cross-site HTML form can't send JSON, so this is what stops forged requests.

Errors are JSON (APIFlask) with a `message` the front end can show.
"""

from typing import Any, Callable, Dict

from apiflask import abort
from flask import Response, jsonify, request
from flask_login import login_required
from sources import models, services
from sources.auxiliary import security_member_workgroup
from sources.extensions import db
from sources.services.inventory import InventoryNotFoundError, InventoryValidationError

from . import inventory_api_bp

EDITOR_ONLY_MESSAGE = (
    "Only principal investigators and senior researchers can change the inventory."
)


def _member_workgroup(workgroup_name: str) -> models.WorkGroup:
    """
    Gets the workgroup if the current user is a member of it.

    Args:
        workgroup_name: The workgroup name from the URL.

    Returns:
        The workgroup.
    """
    workgroup = services.workgroup.from_name(workgroup_name)
    if workgroup is None or not security_member_workgroup(workgroup_name):
        abort(404)
    return workgroup


def _editor_workgroup(workgroup_name: str) -> models.WorkGroup:
    """
    Gets the workgroup if the current user may change its inventory.

    Args:
        workgroup_name: The workgroup name from the URL.

    Returns:
        The workgroup.
    """
    workgroup = _member_workgroup(workgroup_name)
    if not services.inventory.can_edit(workgroup):
        abort(403, EDITOR_ONLY_MESSAGE)
    return workgroup


def _json_body() -> Dict[str, Any]:
    """
    Reads the request's JSON object body. Anything that isn't JSON is refused,
    which is what protects these routes from cross-site form posts.

    Returns:
        The parsed body.
    """
    if not request.is_json:
        abort(415, "Send the request as JSON.")
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        abort(400, "The request body must be a JSON object.")
    return data


def _run(action: Callable[[], Any]) -> Any:
    """
    Runs a service call, turning its errors into JSON error responses and
    discarding any half-made changes.

    Args:
        action: The service call.

    Returns:
        Whatever the service call returns.
    """
    try:
        return action()
    except InventoryValidationError as error:
        db.session.rollback()
        abort(400, str(error))
    except InventoryNotFoundError:
        db.session.rollback()
        abort(404)


# --- Reading -----------------------------------------------------------------


@inventory_api_bp.route("/<workgroup_name>", methods=["GET"])
@login_required
def get_workgroup_inventory(workgroup_name: str) -> Response:
    """
    Gets the workgroup's inventory: its locations, the user's workbooks, the
    user's permissions and every compound it holds with its containers.
    """
    workgroup = _member_workgroup(workgroup_name)
    return jsonify(services.inventory.get_inventory(workgroup))


@inventory_api_bp.route(
    "/<workgroup_name>/items/<int:item_id>/history", methods=["GET"]
)
@login_required
def get_item_history(workgroup_name: str, item_id: int) -> Response:
    """Gets a container's quantity changes, oldest first."""
    workgroup = _member_workgroup(workgroup_name)
    history = _run(lambda: services.inventory.get_item_history(workgroup, item_id))
    return jsonify(history)


@inventory_api_bp.route("/<workgroup_name>/compounds/lookup", methods=["GET"])
@login_required
def lookup_compound(workgroup_name: str) -> Response:
    """
    Looks a compound up by name or CAS number (`?query=`), in the compound
    database and the novel compounds of the user's workbooks.
    """
    workgroup = _member_workgroup(workgroup_name)
    query = request.args.get("query")
    result = _run(lambda: services.inventory.lookup_compound(workgroup, query))
    return jsonify(result)


# --- Containers --------------------------------------------------------------


@inventory_api_bp.route("/<workgroup_name>/items", methods=["POST"])
@login_required
def add_item(workgroup_name: str) -> Response:
    """Adds a container, starting to stock the compound if it's new to the workgroup."""
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    item = _run(lambda: services.inventory.add_item(workgroup, data))
    return jsonify({"id": item.id}), 201


@inventory_api_bp.route(
    "/<workgroup_name>/items/<int:item_id>/adjust", methods=["POST"]
)
@login_required
def adjust_item(workgroup_name: str, item_id: int) -> Response:
    """Sets a container's measured quantity and records the change."""
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    _run(lambda: services.inventory.adjust_item(workgroup, item_id, data))
    return jsonify({"message": "Quantity updated"})


@inventory_api_bp.route("/<workgroup_name>/items/<int:item_id>/move", methods=["POST"])
@login_required
def move_item(workgroup_name: str, item_id: int) -> Response:
    """Moves a container to another location."""
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    _run(lambda: services.inventory.move_item(workgroup, item_id, data))
    return jsonify({"message": "Container moved"})


@inventory_api_bp.route(
    "/<workgroup_name>/items/<int:item_id>/dispose", methods=["POST"]
)
@login_required
def dispose_item(workgroup_name: str, item_id: int) -> Response:
    """Disposes of a container: empties and archives it."""
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    _run(lambda: services.inventory.dispose_item(workgroup, item_id, data))
    return jsonify({"message": "Container disposed"})


# --- Locations ---------------------------------------------------------------


@inventory_api_bp.route("/<workgroup_name>/locations", methods=["POST"])
@login_required
def add_location(workgroup_name: str) -> Response:
    """Adds a storage location."""
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    location = _run(lambda: services.inventory.add_location(workgroup, data))
    return jsonify(services.inventory.serialise_location(location)), 201


@inventory_api_bp.route(
    "/<workgroup_name>/locations/<int:location_id>", methods=["PUT"]
)
@login_required
def edit_location(workgroup_name: str, location_id: int) -> Response:
    """Renames a location or changes its description."""
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    location = _run(
        lambda: services.inventory.edit_location(workgroup, location_id, data)
    )
    return jsonify(services.inventory.serialise_location(location))


@inventory_api_bp.route(
    "/<workgroup_name>/locations/<int:location_id>", methods=["PATCH"]
)
@login_required
def set_location_archived(workgroup_name: str, location_id: int) -> Response:
    """Archives or restores a location (`{"isArchived": true | false}`)."""
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    location = _run(
        lambda: services.inventory.set_location_archived(workgroup, location_id, data)
    )
    return jsonify(services.inventory.serialise_location(location))


# --- Novel compounds ---------------------------------------------------------


@inventory_api_bp.route("/<workgroup_name>/novel-compounds", methods=["POST"])
@login_required
def add_novel_compound(workgroup_name: str) -> Response:
    """
    Adds a novel compound to one of the user's workbooks in this workgroup,
    validated as in the reaction table.
    """
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    compound = _run(lambda: services.inventory.add_novel_compound(workgroup, data))
    return jsonify(compound), 201
