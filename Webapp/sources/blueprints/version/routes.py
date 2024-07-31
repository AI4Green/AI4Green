import contextlib
import os

import toml
from flask import Response, jsonify
from sources import services

from . import version_bp


@version_bp.route("/version", methods=["GET"])
def version() -> Response:
    """
    Returns the current version of the application, with its git commit hash.

    We use the relative path to the root from this file, as there are multiple entrypoints to run the application.
    """
    app_version, git_hash = services.utils.get_app_version()
    return jsonify(f"{app_version}+{git_hash}")
