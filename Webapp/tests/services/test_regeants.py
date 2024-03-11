from flask.testing import FlaskClient
from tests.utils import assert_expected_values, login


def test_reagent_name_non_existent(client: FlaskClient):
    """Tests the response when typing a reagent name that does not exist in the database"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "reagent": "Tes",
        "number": "1",
    }
    # Send a POST request to the route with the form data
    response = client.post("/_reagents", data=form_data)
    assert response.status_code == 200
    expected_values = {
        "reagent": "Tes",
        "reagent_not_found": True,
        "identifiers": [
            "Testamine",
            "Testanol",
            "Testide",
            "Testoic Acid",
            "Testol",
            "Testproductine",
            "Testproductone",
        ],
    }
    assert_expected_values(response.json, expected_values)


def test_reagent_cas_non_existent(client: FlaskClient):
    """Tests the response when typing a reagent CAS that does not exist in the database"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "reagent": "999-09-1",
        "number": "1",
    }
    # Send a POST request to the route with the form data
    response = client.post("/_reagents", data=form_data)
    assert response.status_code == 200
    expected_values = {"reagent": "999-09-1", "cas_not_found": True}
    assert_expected_values(response.json, expected_values)


def test_reagent_compound_table(client: FlaskClient):
    """
    Tests the response when adding a reagent to the reaction table that is in the Compound Table.
    """
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "reagent": "Testide",
        "number": "1",
    }
    # Send a POST request to the route with the form data
    response = client.post("/_reagents", data=form_data)
    # assert status code and response values are as expected
    assert response.status_code == 200
    expected_values = {
        "hazards": "H123",
        "number": "1",
        "concentration": None,
        "density": None,
        "primary_key": 3,
        "name": "Testide",
        "smiles": "CCC",
    }
    assert_expected_values(response.json, expected_values)
