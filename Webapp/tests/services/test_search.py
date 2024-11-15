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
        "mol": """
                  Ketcher 11112414 02D 1   1.00000     0.00000     0

                  1  0  0  0  0  0  0  0  0  0999 V2000
                    6.9750   -3.8500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                M  END""",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "1 results found"
        and response.json.get("status") == "success"
    ), "1 matching result should be found"
    assert (
        response.json.get("images")[0] is not None
    ), "image list should be made"  # TODO: test for actual image generation


def test_structure_search_handler_all_workgroup(client: FlaskClient):
    """Tests the response when looking in all workgroups and books"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "All",
        "workbook": "All",
        "smiles": "C",
        "searchType": "exact_structure",
        "mol": """
                  Ketcher 11112414 02D 1   1.00000     0.00000     0

                  1  0  0  0  0  0  0  0  0  0999 V2000
                    6.9750   -3.8500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                M  END""",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "1 results found"
        and response.json.get("status") == "success"
    ), "1 matching result should be found"
    assert (
        response.json.get("images")[0] is not None
    ), "image list should be made"  # TODO: test for actual image generation


def test_structure_search_handler_all_workbooks(client: FlaskClient):
    """Tests the response when looking in all workbooks for a specific workgroup"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "All",
        "smiles": "C",
        "searchType": "exact_structure",
        "mol": """
                  Ketcher 11112414 02D 1   1.00000     0.00000     0

                  1  0  0  0  0  0  0  0  0  0999 V2000
                    6.9750   -3.8500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                M  END""",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "1 results found"
        and response.json.get("status") == "success"
    ), "1 matching result should be found"
    assert (
        response.json.get("images")[0] is not None
    ), "image list should be made"  # TODO: test for actual image generation


def test_structure_search_handler_no_results(client: FlaskClient):
    """Tests the response when there are no matches found in the search"""

    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "smiles": "C1N2CN3CN1CN(C2)C3",
        "searchType": "exact_structure",
        "mol": """
              Ketcher 11112414212D 1   1.00000     0.00000     0

             10 12  0  0  0  0  0  0  0  0999 V2000
                9.7566   -3.8684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               10.2893   -4.6405    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
               11.2243   -4.5652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               11.6266   -3.7178    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
               11.0939   -2.9457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               10.1589   -3.0210    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
                9.2906   -3.1441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                8.9234   -3.9789    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
                9.5005   -4.7043    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               10.8924   -3.3821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
              1  2  1  0     0  0
              2  3  1  0     0  0
              3  4  1  0     0  0
              4  5  1  0     0  0
              5  6  1  0     0  0
              6  1  1  0     0  0
              6  7  1  0     0  0
              7  8  1  0     0  0
              8  9  1  0     0  0
              9  2  1  0     0  0
              8 10  1  0     0  0
             10  4  1  0     0  0
            M  END
            """,
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "No results found"
        and response.json.get("status") == "fail"
    ), "No matching results should be found"
    assert response.json.get("images") is None, "image list should not be made"


def test_structure_search_handler_invalid_smiles(client: FlaskClient):
    """Tests the response when there are no matches found in the search"""

    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "smiles": "thisIsDefinitelyNotASmilesString",
        "searchType": "exact_structure",
        "mol": "",
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "Invalid structure, please try again!"
        and response.json.get("status") == "fail"
    ), "No matching results should be found"
    assert response.json.get("images") is None, "image list should not be made"


def test_structure_search_handler_polymer(client: FlaskClient):
    """Tests the response when there is a polymer found in the search"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "smiles": "*C*",
        "searchType": "exact_structure",
        "mol": """
                  Ketcher 11142416132D 1   1.00000     0.00000     0

                  3  2  0  0  0  0  0  0  0  0999 V2000
                    9.6721    4.5103    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
                   10.6721    4.5103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                   11.1721    3.6443    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
                  1  2  1  0     0  0
                  2  3  1  0     0  0
                M  STY  1   1 SRU
                M  SLB  1   1   1
                M  SCN  1   1 HT
                M  SMT   1 n
                M  SAL   1  1   2
                M  SBL   1  2   1   2
                M  SDI   1  4   10.1721    4.0103   10.1721    5.0103
                M  SDI   1  4   11.3551    4.3273   10.4890    3.8273
                M  SDS EXP  1   1
                M  END
                """,
    }
    # Send a POST request to the route with the form data
    response = client.post("/structure_search_handler", data=form_data)
    # Assert that the response is as expected
    assert response.status_code == 200
    assert (
        response.json.get("message") == "1 results found"
        and response.json.get("status") == "success"
    ), "1 matching result should be found"
    assert response.json.get("images")[0] is not None, "image list should be made"
