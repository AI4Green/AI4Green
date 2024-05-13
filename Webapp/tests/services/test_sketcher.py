from flask.testing import FlaskClient
from tests.utils import login


def test_sketcher(client: FlaskClient):
    """Tests we can load the sketcher page"""
    login(client)
    # make a request and assert it is successful and as expected
    url = make_url()
    response = client.post(url)
    assert response.status_code == 200
    # assert sketchers are present
    assert (
        b"""<div id="ketcher-sketcher"></div>\n        <div id="marvin-sketcher" style="width: 1020px; height: 480px"></div>"""
        in response.data
    )


def test_demo_sketcher(client: FlaskClient):
    """Tests we can load the sketcher page"""
    login(client)
    # make a request and assert it is successful and as expected
    response = client.post("demo")
    assert response.status_code == 200
    # assert sketchers are present
    assert (
        b"""<div id="ketcher-sketcher"></div>\n        <div id="marvin-sketcher" style="width: 1020px; height: 480px"></div>"""
        in response.data
    )


def make_url(
    workgroup: str = "Test-Workgroup",
    workbook: str = "Test-Workbook",
    reaction_id: str = "TW1-001",
    tutorial: str = "no",
) -> str:
    """
    Constructs a URL for making a GET request based on tutorial, workgroup, workbook, and reaction ID.

    Args:
        workgroup (str, optional): The workgroup associated with the reaction. Defaults to "Test-Workgroup".
        workbook (str, optional): The workbook associated with the reaction. Defaults to "Test-Workbook".
        reaction_id (str, optional): The ID of the reaction. Defaults to "TW1-001".
        tutorial (str, optional): Whether the reaction is being performed in tutorial mode. Defaults to "no".

    Returns:
        str: The constructed URL for making a GET request.
    """
    # Defining the endpoint
    endpoint = "/sketcher"
    # Constructing the URL with the arguments
    url = f"{endpoint}/{workgroup}/{workbook}/{reaction_id}/{tutorial}"
    return url
