import json
from typing import Dict, Tuple

import requests
from requests import Response


def retrosynthesis_api_call(request_url: str, retrosynthesis_base_url: str) -> str:
    """Make a call to the retrosynthesis API to trigger a retrosynthesis job.

    Args:
        request_url (str): the url containing the retrosynthesis server and target SMILES
        retrosynthesis_base_url (str): the base url for the retrosynthesis server

    Returns:
        Tuple[str, str, Union[str, List[Dict]]]: The job ID which can be polled for results
    """
    try:
        response = requests.get(request_url)
    except requests.exceptions.ConnectionError:
        return "Retrosynthesis server is down. If this problem persists please report this error to admin@ai4green.app"
    if response.status_code == 500:
        return "Error with molecule"

    try:
        response_data = response.json()
        job_id = response_data["job_id"]
        return job_id
    except Exception:
        failure_message = determine_retrosynthesis_error(retrosynthesis_base_url)
        return "failed", failure_message, ""


def retrosynthesis_results_poll(request_url: str) -> Tuple[str, dict, dict]:
    """Poll the retrosynthesis API for the results of a job.

    Args:
        request_url (str): The URL to the results endpoint, containing the job ID.

    Returns:
        Tuple[str, dict, dict]: Status, solved routes, raw routes
    """
    try:
        response = requests.get(request_url)
        if not response.ok:
            return "error", {}, {}

        data = response.json()
        status = data.get("status", "error")
        if status == "error":
            return status, {}, {}

        results = data.get("results", {})
        solved_routes = results.get("solved_route_dict", {})
        raw_routes = results.get("raw_routes", {})
        return status, solved_routes, raw_routes
    except Exception:
        return "error", {}, {}


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
