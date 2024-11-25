from flask import Flask
from flask.testing import FlaskClient
from sources import services
from tests.utils import assert_expected_values, login


def test_novel_compound_sketcher(app: Flask, client: FlaskClient):
    """Tests the response when adding a novel compound from the sketcher"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "name": "A Novel compound",
        "source": "sketcher",
        "smiles": "CCN",
        "density": "",
        "concentration": "",
        "hPhrase": "",
        "cas": "",
        "molWeight": "45.08",
        "component": "component",
    }
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_novel_compound", data=form_data)
    assert (
        response.status_code == 200
        and response.json["feedback"] == "Compound added to the database"
    )
    # confirm compound is in database
    with app.app_context():
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            "Test-Workgroup", "Test-Workbook"
        )

        new_compound = services.novel_compound.from_inchi_and_workbook(
            "InChI=1S/C2H7N/c1-2-3/h2-3H2,1H3", workbook
        )
        assert new_compound.name == "A Novel compound"


def test_novel_compound_solvent(app: Flask, client: FlaskClient):
    """Tests the response when adding a novel compound from the sketcher"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "name": "A Novel solvent",
        "source": "table",
        "smiles": "",
        "density": "",
        "concentration": "",
        "hPhrase": "",
        "cas": "",
        "molWeight": "",
        "component": "solvent",
    }
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_novel_compound", data=form_data)
    assert (
        response.status_code == 200
        and response.json["feedback"] == "Compound added to the database"
    )
    # confirm compound is in database, has the right flag, and is linked to a novel compound table entry
    with app.app_context():
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            "Test-Workgroup", "Test-Workbook"
        )

        solvent_list = services.solvent.get_workbook_list(workbook)
        novel_solvent = [
            solvent for solvent in solvent_list if solvent.name == "A Novel solvent"
        ][0]
        assert novel_solvent.flag == 5
        assert novel_solvent.novel_compound is not None


def test_novel_compound_table_all_data(app: Flask, client: FlaskClient):
    """Tests the response when adding a novel compound from the table with all the data"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "name": "A Novel Reagent",
        "source": "table",
        "smiles": "CCO",
        "density": "0.95",
        "concentration": "",
        "hPhrase": "H900-H901",
        "cas": "874-43-1",
        "molWeight": "46.068",
        "component": "component",
    }
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_novel_compound", data=form_data)
    assert (
        response.status_code == 200
        and response.json["feedback"] == "Compound added to the database"
    )
    # confirm compound is in database
    with app.app_context():
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            "Test-Workgroup", "Test-Workbook"
        )

        new_compound = services.novel_compound.from_inchi_and_workbook(
            "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", workbook
        )
        assert new_compound.name == "A Novel Reagent"
        assert new_compound.hphrase == "H900-H901"
        assert new_compound.cas == "874-43-1"
        assert new_compound.density == 0.95
