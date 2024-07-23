import json
from typing import Dict, List, Tuple, Union

import requests
from requests import Response


def retrosynthesis_api_call(
    request_url: str, retrosynthesis_base_url: str
) -> Tuple[str, str, Union[str, List[Dict]]]:
    """
    Makes a call to the retrosynthesis api server with the target molecule as a SMILES strings and recieves a JSON
    containing the retrosynthetic route
    Args:
        request_url - the url containing the retrosynthesis server and target SMILES
        retrosynthesis_base_url - the base url for the retrosynthesis server
    Returns:
         a string to indicate if the api call failed or succeeded
         the user feedback message either with the error type or a success message
         The list of routes as dictionaries
    """

    try:
        response = requests.get(request_url)
    except requests.exceptions.ConnectionError:
        return (
            "failed",
            "Retrosynthesis server is down. If this problem persists please report this error to admin@ai4green.app",
            "",
        )
    if response.status_code == 500:
        return "failed", "Error with molecule", ""

    try:
        solved_routes = json.loads(response.content)["Message"]
        if solved_routes == "invalid key":
            return "failed", "invalid key", ""
        if solved_routes == {}:
            return (
                "failed",
                "Could not find a successful route to the molecule.",
                "",
            )
        assert solved_routes
        # if query fails try and print error code but
    except Exception:
        failure_message = determine_retrosynthesis_error(retrosynthesis_base_url)
        return "failed", failure_message, ""
    validation, validation_message = validate_retrosynthesis_api_response(
        response, solved_routes
    )
    return validation, validation_message, solved_routes


def validate_retrosynthesis_api_response(
    response: Response, solved_routes: Dict
) -> Tuple[str, str]:
    """
    Assesses whether the api response was successful or a failure
    Args:
        response - the response from the retrosynthesis api server
        solved_routes - the retrosynthetic routes from the server
    Returns:
        string to indicate failed or success of the api call
        message describing why route was invalid if that was the case

    """
    if response.status_code == 200 and solved_routes:
        return "success", ""
    if not solved_routes:
        return "failed", "No routes found for molecule"
    if response:
        if response.status_code != 200:
            return (
                "failed",
                "Error with retrosynthesis function. Please try again in a few minutes. "
                "If this is a repeated error please report this to admin@ai4green.app, with details of your error",
            )


def determine_retrosynthesis_error(retrosynthesis_base_url: str) -> str:
    """
    Tests the server to try and determine the error
    Args:
        retrosynthesis_base_url - the base url of the retrosynthesis server
    Returns:
        a message describing the error
    """
    url = f"{retrosynthesis_base_url}/"
    response_content = ""
    try:
        response = requests.get(url)
        response_content = json.loads(response.content["Message"])
    except Exception:
        return "Retrosynthesis server is down. If this problem persists please report this error to admin@ai4green.app"
    if not response_content == "Retrosynthesis service is running":
        return (
            "retrosynthesis service is currently having technical issues. If this problem persists please report this "
            "error to admin@ai4green.app"
        )
    return "Error with retrosynthesis. Please report this error to admin@ai4green.app"
