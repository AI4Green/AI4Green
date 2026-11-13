"""
Search API for the React search page.

Access:
- Every route needs a logged-in user, and only searches the workbooks that
  user belongs to. Asking for a workgroup or workbook they can't see gets a
  404, so the routes don't reveal which ones exist.
- The search is a POST because a drawn structure can be too long for a URL.
  It only reads data, but it still only accepts a JSON body, in line with the
  rest of the API's protection against cross-site form posts.

Errors are JSON (APIFlask) with a `message` the front end can show.
"""

from typing import Any, Dict

from apiflask import abort
from flask import Response, jsonify, request
from flask_login import login_required
from sources import services
from sources.services.reaction_search import (
    ReactionSearchScopeError,
    ReactionSearchValidationError,
)
from sources.services.structure_search import StructureSearchError

from . import search_api_bp


def _json_body() -> Dict[str, Any]:
    """
    Reads the request's JSON object body, refusing anything that isn't JSON.

    Returns:
        The parsed body.
    """
    if not request.is_json:
        abort(415, "Send the request as JSON.")
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        abort(400, "The request body must be a JSON object.")
    return data


@search_api_bp.route("/workgroups", methods=["GET"])
@login_required
def get_search_workgroups() -> Response:
    """Gets the user's workgroups and their workbooks, for the search filters."""
    return jsonify({"workgroups": services.reaction_search.accessible_workbooks()})


@search_api_bp.route("/reactions", methods=["POST"])
@login_required
def search_reactions() -> Response:
    """
    Searches the user's reactions. The body may contain `workgroup`,
    `workbook`, `keyword`, `reagent`, `dateFrom`, `dateTo` (YYYY-MM-DD) and
    `structure` (`{smiles, mode: "exact" | "substructure", roles}`), and needs
    at least one of the search terms.

    Returns `{results: [...], truncated}`.
    """
    data = _json_body()
    try:
        criteria = services.reaction_search.read_criteria(data)
        results = services.reaction_search.search_reactions(criteria)
    except (ReactionSearchValidationError, StructureSearchError) as error:
        abort(400, str(error))
    except ReactionSearchScopeError:
        abort(404)
    return jsonify(services.reaction_search.serialise_results(results))
