import re
from typing import Dict


def update_dict_keys_to_camelCase(snake_dict: Dict[str, str]) -> Dict[str, str]:
    camel_case_dict = {
        snake_case_to_camelCase(key): value for key, value in snake_dict.items()
    }
    return camel_case_dict


def update_dict_keys_to_snake_case(camelCase_dict: Dict[str, str]) -> Dict[str, str]:
    snake_case_dict = {
        camelCase_to_snake_case(key): value for key, value in camelCase_dict.items()
    }
    return snake_case_dict


def snake_case_to_camelCase(string: str) -> str:
    """
    Converts a string from snake_case to camelCase which is often used in JavaScript.
    """
    parts = re.split(r"[ _]+", string)  # Split by space or underscore
    # Capitalize the first character of each part except the first one
    camel_case = "".join(
        [part.capitalize() if i > 0 else part for i, part in enumerate(parts)]
    )
    return camel_case


def camelCase_to_snake_case(string: str) -> str:
    """
    Converts a string from camelCase which is often used in JavaScript to the preferred snake_case.
    """
    return "".join(["_" + c.lower() if c.isupper() else c for c in string]).lstrip("_")
