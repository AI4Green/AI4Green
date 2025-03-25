from flask import Flask, json
from flask.testing import FlaskClient
from sources import services
from sources.services.polymer_novel_compound import find_canonical_repeats
from tests.utils import login


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


def test_polymer_novel_compound_sketcher(app: Flask, client: FlaskClient):
    """Tests the response when adding a novel compound from the sketcher"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "name": "A Polymer Novel compound",
        "source": "sketcher",
        "smiles": json.dumps(["PC(CC(C)P)C", "C"]),
        "density": "",
        "concentration": "",
        "hPhrase": "",
        "cas": "",
        "molWeight": json.dumps(["45.08", "10"]),
        "component": "component",
    }
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_polymer_novel_compound", data=form_data)
    assert (
        response.status_code == 200
        and response.json["feedback"] == "Compound added to the database"
    )
    # confirm compound is in database
    with app.app_context():
        workbook = services.workbook.get_workbook_from_group_book_name_combination(
            "Test-Workgroup", "Test-Workbook"
        )

        new_compound = services.polymer_novel_compound.from_smiles_and_workbook(
            ["PC(CC(C)P)C", "C"], workbook
        )
        assert new_compound.name == "A Polymer Novel compound"


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


def test_polymer_smiles_parsing():
    """Tests if different polymers are parsed correctly"""
    assert find_canonical_repeats("CC{-}C{+n}C") == ["*C*"]  # end groups

    # different branching patterns
    assert find_canonical_repeats("*C{-}(CC{+n}*)C") == ["*CCC(*)C"]
    assert find_canonical_repeats("*C{-}C{+n}(C)(C)*") == ["*CC(*)(C)C"]
    assert find_canonical_repeats("C(C)(C)C{-}(CC{+n}C(C)C)C") == ["*CCC(*)C"]
    assert find_canonical_repeats("*C{-}(C(C(C{+n}*)=O)C(C)C)=O") == [
        "*C(=O)CC(=O)C(*)C(C)C"
    ]
    assert find_canonical_repeats("*C{-}C(C(C{+n}*)=O)C(C)C") == ["*CCC(=O)C(*)C(C)C"]
    assert find_canonical_repeats(
        "*C1{-}C(C23C1C2(C3C{+n}(C(C)Cl)C(C)C)C(Cl)C)C1CC1"
    ) == ["*C(C(C)Cl)C1C2(C(C)Cl)C3C(*)C(C4CC4)C132"]

    assert find_canonical_repeats("*C1{-}CC{+n}(CCC1)*") == ["*C1CCCC(*)C1"]  # rings
    assert find_canonical_repeats("C[SiH2]{-}CC{+n}C") == [
        "*C[SiH2]C*"
    ]  # groups in [] brackets
    assert find_canonical_repeats("C(C{-}{+n}*)C{-}(C{+n}*)C") == [
        "*C*",
        "*CC(*)C",
    ]  # copolymers

    assert find_canonical_repeats("C[*]{-}C{+n}C") == "dummy"  # dummy atoms = blocked
