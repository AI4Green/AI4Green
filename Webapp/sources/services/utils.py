import contextlib
import datetime
from audioop import error
from typing import Dict, Tuple, Union

import flask
import ipinfo
import toml
from flask import current_app, request


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


def get_privacy_policy_date() -> datetime.datetime:
    """
    Returns the date the current privacy policy was introduced in to the app.

    Returns:
        datetime.datetime, datetime the privacy policy was introduced in to the app.
    """
    return datetime.datetime(2025, 3, 18)


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


def get_ip_address() -> Union[str, None]:
    """Returns current IP address"""
    if request.headers.get("X-Forwarded-For"):
        ip = request.headers.get("X-Forwarded-For").split(",")[0].strip()
        return ip.split(":")[0]
    return None


def get_location() -> Dict[str, str]:
    """
    Uses IPInfo call to get user location

    Returns:
        Dict[str, str]: dictionary with IP address, country name and city of the request
    """
    handler = ipinfo.getHandler(current_app.config["IPINFO_API_KEY"])
    ip = get_ip_address()
    if ip is not None:
        location = handler.getDetails(ip)
        return {
            "IP_address": ip,
            "country": location.country_name,
            "city": location.city,
        }
    # location information is not supported for local instances
    return {
        "IP_address": "Location Not Found",
        "country": "Location Not Found",
        "city": "Location Not Found",
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
