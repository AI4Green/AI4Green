import contextlib
from typing import Tuple, Dict

import toml
from flask import current_app, request
import ipinfo


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


def check_positive_number(s: float) -> bool:
    """Checks the entry is a positive number."""
    try:
        return s >= 0
    except (ValueError, TypeError):
        return False


def remove_duplicates_keep_first(lst: list) -> list:
    """
    Removes duplicates from a list but keeps the order and only removes the non-first duplicates
    Args:
        lst - a list potentially containing duplicates
    Returns:
        a list without duplicates
    """
    seen_items = set()
    new_list = []

    for item in lst:
        if item not in seen_items:
            new_list.append(item)
            seen_items.add(item)

    return new_list


def get_ip_address() -> str:
    """Returns current IP address"""
    return request.remote_addr


def get_location() -> Dict[str, str]:
    handler = ipinfo.getHandler(current_app.config["IPINFO_API_KEY"])
    ip = get_ip_address()
    location = handler.getDetails(ip)
    try:
        # This will fail if running on local host
        return {
            "IP_address": ip,
            "country": location.country_name,
            "city": location.city,
        }
    except AttributeError:
        # If running on local host use default values that point to Nottingham
        print("Local Host IP is not findable with IPInfo. Using Nottingham as default...")
        return {
            "IP_address": "NOTTINGHAM_IP_ADDRESS",
            "country": "United Kingdom",
            "city": "Nottingham",
        }


def remove_spaces_and_dashes(name: str) -> str:
    """
    Removes spaces and dashes from input string
    Args:
        name: str, the name of the group which needs white space removed.
    Returns:
        name: str, the name without spaces.
    """
    return name.replace(" ", "").replace("-", "")
