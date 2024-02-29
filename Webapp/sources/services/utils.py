def camelCase_to_snake_case(string: str) -> str:
    """
    Converts a string from camelCase which is often used in JavaScript to the preferred snake_case.
    """
    return "".join(["_" + c.lower() if c.isupper() else c for c in string]).lstrip("_")
