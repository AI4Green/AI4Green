"""
Inventory API for the React inventory pages.

Access:
- Every route needs a logged-in member of the workgroup. Anyone else gets a
  404, so the routes don't reveal which workgroups exist.
- Reading is open to all members; changes need a principal investigator or
  senior researcher (services.inventory.EDITOR_ROLES), otherwise 403.
- Writes only accept JSON bodies. App-wide CSRF tokens are disabled, and a
  cross-site HTML form can't send JSON, so this is what stops forged requests.
  The CSV import has to accept a multipart file upload, which a form *can*
  send, so it requires the IMPORT_HEADER custom header instead: a cross-site
  page can't add one without a CORS preflight, which this API doesn't allow.

Errors are JSON (APIFlask) with a `message` the front end can show.
"""

from typing import Any, Callable, Dict

from apiflask import abort
from flask import Response, jsonify, request
from flask_login import current_user, login_required
from sources import models, services
from sources.extensions import db
from sources.services.cheminventory import (
    ChemInventoryError,
    ChemInventoryNotConnectedError,
    ChemInventoryRejectedError,
    ChemInventoryUnavailableError,
    ChemInventoryValidationError,
)
from sources.services.inventory import InventoryNotFoundError, InventoryValidationError

from . import inventory_api_bp

EDITOR_ONLY_MESSAGE = (
    "Only principal investigators and senior researchers can change the inventory."
)
# Custom header the CSV import routes require, as their CSRF protection
IMPORT_HEADER = "X-AI4Green-Request"
IMPORT_HEADER_VALUE = "inventory-import"


def _member_workgroup(workgroup_name: str) -> models.WorkGroup:
    """
    Gets the workgroup if the current user is a member of it.

    The decorators in sources.decorators make the same check, but they flash a
    message and redirect to the home page, which a JSON API can't use.

    Args:
        workgroup_name: The workgroup name from the URL.

    Returns:
        The workgroup.
    """
    workgroup = services.workgroup.from_name(workgroup_name)
    if (
        workgroup is None
        or services.workgroup.get_user_type(workgroup_name, current_user) is None
    ):
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


def _require_import_header() -> None:
    """
    Refuses multipart requests without the import's custom header, which is
    what protects the CSV import from cross-site form posts.
    """
    if request.headers.get(IMPORT_HEADER) != IMPORT_HEADER_VALUE:
        abort(400, "This request is missing a required header.")


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


def _run_cheminventory(action: Callable[[], Any], token_check: bool = False) -> Any:
    """
    Runs a ChemInventory service call, turning its errors into JSON error
    responses. Every message comes from the service and is safe to show.

    Args:
        action: The service call.
        token_check: Whether the call is checking a token the user just
            entered, in which case ChemInventory refusing it is the user's
            input error (400) rather than a problem with ChemInventory (502).

    Returns:
        Whatever the service call returns.
    """
    try:
        return action()
    except ChemInventoryValidationError as error:
        db.session.rollback()
        abort(400, str(error))
    except ChemInventoryRejectedError as error:
        db.session.rollback()
        abort(400 if token_check else 502, str(error))
    except ChemInventoryError as error:
        db.session.rollback()
        abort(502, str(error))
    except ChemInventoryNotConnectedError as error:
        db.session.rollback()
        abort(409, str(error))
    except ChemInventoryUnavailableError as error:
        db.session.rollback()
        abort(503, str(error))


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


@inventory_api_bp.route("/<workgroup_name>/structure-search", methods=["POST"])
@login_required
def structure_search(workgroup_name: str) -> Response:
    """
    Finds the workgroup's compounds matching a drawn structure. The body is
    `{smiles, mode}`, where mode is "exact" or "substructure". It's a POST
    because structures can be too long for a URL; it only reads data.
    """
    workgroup = _member_workgroup(workgroup_name)
    data = _json_body()
    ids = _run(lambda: services.inventory.search_by_structure(workgroup, data))
    return jsonify({"workgroupCompoundIds": ids})


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


# --- CSV import --------------------------------------------------------------
# Both routes take a multipart upload (`file`, optional `workbookId`) and run
# the whole import. Preview always rolls back; import saves only if every row
# is valid. Neither keeps anything between requests.


def _import_csv(workgroup_name: str, save: bool) -> Dict[str, Any]:
    """
    Checks access and the CSRF header, then runs the CSV import.

    Args:
        workgroup_name: The workgroup name from the URL.
        save: Whether to save the import, or only preview it.

    Returns:
        The import result from services.inventory_import.import_csv.
    """
    workgroup = _editor_workgroup(workgroup_name)
    _require_import_header()
    file = request.files.get("file")
    workbook_id = request.form.get("workbookId")
    return _run(
        lambda: services.inventory_import.import_csv(
            workgroup, file, workbook_id, save=save
        )
    )


@inventory_api_bp.route("/<workgroup_name>/import/preview", methods=["POST"])
@login_required
def preview_import(workgroup_name: str) -> Response:
    """Checks a CSV import and shows what it would do, without saving it."""
    return jsonify(_import_csv(workgroup_name, save=False))


@inventory_api_bp.route("/<workgroup_name>/import", methods=["POST"])
@login_required
def import_inventory(workgroup_name: str) -> Response:
    """
    Imports a CSV file of containers. Saves nothing unless every row is valid;
    otherwise returns 400 with the same per-row results as the preview.
    """
    result = _import_csv(workgroup_name, save=True)
    status = 201 if result["saved"] else 400
    return jsonify(result), status


# --- ChemInventory -----------------------------------------------------------
# Read-only view of a workgroup's ChemInventory account. Any member can view
# it; principal investigators and senior researchers can set or remove the
# API token, which is never sent back to the browser.


@inventory_api_bp.route("/<workgroup_name>/cheminventory/search", methods=["POST"])
@login_required
def search_cheminventory(workgroup_name: str) -> Response:
    """
    Searches the workgroup's ChemInventory by CAS number or name
    (`{"type": "cas" | "name", "query": ...}`), returning the matching
    containers grouped by location. A JSON POST rather than a GET, so another
    site can't make a member's browser run searches (e.g. with an image tag)
    to load ChemInventory and this server.
    """
    workgroup = _member_workgroup(workgroup_name)
    data = _json_body()
    results = _run_cheminventory(lambda: services.cheminventory.search(workgroup, data))
    return jsonify(results)


@inventory_api_bp.route("/<workgroup_name>/cheminventory/connection", methods=["GET"])
@login_required
def get_cheminventory_connection(workgroup_name: str) -> Response:
    """Gets whether the workgroup has connected ChemInventory, without the token."""
    workgroup = _member_workgroup(workgroup_name)
    return jsonify(services.cheminventory.get_connection_status(workgroup))


@inventory_api_bp.route("/<workgroup_name>/cheminventory/connection", methods=["PUT"])
@login_required
def save_cheminventory_connection(workgroup_name: str) -> Response:
    """
    Connects ChemInventory with an API token (`{"token": ...}`), or replaces
    the token. The token is checked with ChemInventory before it's saved.
    """
    workgroup = _editor_workgroup(workgroup_name)
    data = _json_body()
    _run_cheminventory(
        lambda: services.cheminventory.save_connection(workgroup, data),
        token_check=True,
    )
    return jsonify(services.cheminventory.get_connection_status(workgroup))


@inventory_api_bp.route(
    "/<workgroup_name>/cheminventory/connection", methods=["DELETE"]
)
@login_required
def delete_cheminventory_connection(workgroup_name: str) -> Response:
    """
    Disconnects ChemInventory, deleting the stored token. Takes no body: a
    cross-site DELETE needs a CORS preflight, which this app never allows, so
    another site can't send it.
    """
    workgroup = _editor_workgroup(workgroup_name)
    _run_cheminventory(lambda: services.cheminventory.delete_connection(workgroup))
    return jsonify({"message": "ChemInventory disconnected"})
