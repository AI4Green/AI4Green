import contextlib
import os

import toml
from flask import Response, jsonify

from . import version_bp


@version_bp.route("/version", methods=["GET"])
def version() -> Response:
    """
    Returns the current version of the application, with its git commit hash.

    We use the relative path to the root from this file, as there are multiple entrypoints to run the application.
    """
    app_directory = os.path.dirname(os.path.abspath(__file__))

    # Traverse up directories to reach the project root
    project_root = os.path.abspath(os.path.join(app_directory, "..", "..", "..", ".."))

    # Construct paths relative to the project root
    pyproject_path = os.path.join(project_root, "pyproject.toml")
    hash_file_path = os.path.join(project_root, "hash.txt")

    app_version = ""
    with contextlib.suppress(FileNotFoundError):
        with open(pyproject_path, "r") as f:
            config = toml.load(f)
            app_version = config["tool"]["poetry"]["version"]

    git_hash = ""
    with contextlib.suppress(FileNotFoundError):
        with open(hash_file_path, "r") as f:
            git_hash = f.read().strip()

    return jsonify(f"{app_version}+{git_hash}")
