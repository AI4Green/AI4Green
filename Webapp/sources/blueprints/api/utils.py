"""Helpers shared by the JSON API blueprints."""

from typing import Any, Dict, List, Union

from apiflask import abort
from flask import request


def json_body(allow_list: bool = False) -> Union[Dict[str, Any], List[Any]]:
    """
    Reads the request's JSON body. Anything that isn't JSON is refused, which is
    what protects write routes from cross-site form posts (app-wide CSRF tokens are
    disabled for the API).

    Args:
        allow_list: Accept a JSON array as well as an object.

    Returns:
        The parsed body.
    """
    if not request.is_json:
        abort(415, "Send the request as JSON.")
    data = request.get_json(silent=True)
    allowed_types = (dict, list) if allow_list else (dict,)
    if not isinstance(data, allowed_types):
        abort(400, "The request body must be a JSON object.")
    return data
