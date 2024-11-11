import os

from flask.testing import FlaskClient
from tests.utils import login

IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


def test_structure_search_handler_success(client: FlaskClient):
    """Tests the response when there is a match found in the search"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "smiles": "C",
        "searchType": "exact_structure",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "1 results found"
        and response.json.get("status") == "success"
    ), "1 matching result should be found"
    assert "svg" in response.json.get("schemes")[0], "SVG image should be made"


def test_structure_search_handler_all_workgroup(client: FlaskClient):
    """Tests the response when looking in all workgroups and books"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "All",
        "workbook": "All",
        "smiles": "C",
        "searchType": "exact_structure",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "1 results found"
        and response.json.get("status") == "success"
    ), "1 matching result should be found"
    assert "svg" in response.json.get("schemes")[0], "SVG image should be made"


def test_structure_search_handler_all_workbooks(client: FlaskClient):
    """Tests the response when looking in all workbooks for a specific workgroup"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "All",
        "smiles": "C",
        "searchType": "exact_structure",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "1 results found"
        and response.json.get("status") == "success"
    ), "1 matching result should be found"
    assert "svg" in response.json.get("schemes")[0], "SVG image should be made"


def test_structure_search_handler_no_results(client: FlaskClient):
    """Tests the response when there are no matches found in the search"""

    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "smiles": "C1N2CN3CN1CN(C2)C3",
        "searchType": "exact_structure",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "No results found"
        and response.json.get("status") == "fail"
    ), "No matching results should be found"
    assert response.json.get("schemes") is None, "SVG image should not be made"


def test_structure_search_handler_invalid_smiles(client: FlaskClient):
    """Tests the response when there are no matches found in the search"""

    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "smiles": "thisIsDefinitelyNotASmilesString",
        "searchType": "exact_structure",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "Invalid structure, please try again!"
        and response.json.get("status") == "fail"
    ), "No matching results should be found"
    assert response.json.get("schemes") is None, "SVG image should not be made"
