from flask.testing import FlaskClient
from tests.utils import assert_expected_values, login


def test_solvent_non_existent(client: FlaskClient):
    """Tests the response when typing a solvent name that does not exist in the database"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "solvent": "Tea",
        "number": "1",
    }
    # Send a POST request to the route with the form data
    response = client.post("/_solvents", data=form_data)
    assert response.status_code == 204


def test_solvent_chem21(client: FlaskClient):
    """
    Tests the response when adding a solvent to the reaction table that has a CHEM21 flag and entry in the solvent table
    """
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "solvent": "Testanol",
        "number": "1",
    }
    # Send a POST request to the route with the form data
    response = client.post("/_solvents", data=form_data)
    # assert status code and response values are as expected
    assert response.status_code == 200
    expected_values = {
        "alert_message": "",
        "flag": "hazard-hazardous",
        "hazards": "H900",
        "num": "1",
        "primary_key": [5],
        "new_solvent": False,
        "solvent": "Testanol",
    }
    assert_expected_values(response.json, expected_values)
