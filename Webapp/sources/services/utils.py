import contextlib
import os
import re
from typing import Tuple

import toml


def camelCase_to_snake_case(string: str) -> str:
    """
    Converts a string from camelCase which is often used in JavaScript to the preferred snake_case.
    """
    return "".join(["_" + c.lower() if c.isupper() else c for c in string]).lstrip("_")


def get_app_version() -> Tuple[str, str]:
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

    return app_version, git_hash


def check_is_number(string: str) -> bool:
    """Checks whether a string is a numerical value including decimals"""
    try:
        float(string)
        return True
    except ValueError:
        return False
