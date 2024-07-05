import contextlib
from typing import Tuple

import toml
from flask import current_app


def camelCase_to_snake_case(string: str) -> str:
    """
    Converts a string from camelCase which is often used in JavaScript to the preferred snake_case.
    """
    return "".join(["_" + c.lower() if c.isupper() else c for c in string]).lstrip("_")


def get_app_version() -> Tuple[str, str]:
    """
    Returns the current version of the application, with its git commit hash.

    We use the relative path to the root from this file, as there are multiple entrypoints to run the application.

    Returns:
        app_version - the current version of the app as specified in the pyproject.toml
        git_hash - the git commit hash saved in hash.txt
    """
    app_version = ""
    with contextlib.suppress(FileNotFoundError):
        with open(current_app.config["PYPROJECT_PATH"], "r") as f:
            config = toml.load(f)
            app_version = config["tool"]["poetry"]["version"]

    git_hash = ""
    with contextlib.suppress(FileNotFoundError):
        with open(current_app.config["HASH_FILE_PATH"], "r") as f:
            git_hash = f.read().strip()

    return app_version, git_hash


def check_is_number(string: str) -> bool:
    """Checks whether a string is a numerical value including decimals"""
    try:
        float(string)
        return True
    except ValueError:
        return False
