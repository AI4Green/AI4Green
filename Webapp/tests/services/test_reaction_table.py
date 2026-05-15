import json
from typing import List, Optional

from bs4 import BeautifulSoup
from flask.testing import FlaskClient
from tests.utils import assert_expected_values, login


def test_autoupdate_reaction_table_compounds_present(client: FlaskClient):
    login(client)

    # compounds are in the db
    data = get_test_data(reaction_smiles="C.CC>>CP.CCP")

    response = client.post("/autoupdate_reaction_table", json=data)

    assert_reaction_table_response_for_test_compounds(response)


def test_autoupdate_reaction_table_ionic_reactions(client: FlaskClient):
    login(client)

    data = get_test_data(
        reaction_smiles="CC(C(=O)N1CCCC1C(=O)OC(C)(C)C)[N+]([O-])=O.F[B-](F)(F)F.FC1=CC=C(C=C1)[Bi+](C1=CC=C(F)C=C1)(C1=CC=C(F)C=C1)C1=CC=C(F)C=C1>>CC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[N+]([O-])=O |f:1.2|"
    )

    response = client.post("/autoupdate_reaction_table", json=data)

    # ensure that novel compound template is rendered
    assert response.status_code == 200
    assert (
        response.json["novelCompound"] is True
        and "Reactant 1 not in database" in response.json["reactionTable"]
    )

    # replace first compound with a db compound
    data2 = get_test_data(
        reaction_smiles="C.F[B-](F)(F)F.FC1=CC=C(C=C1)[Bi+](C1=CC=C(F)C=C1)(C1=CC=C(F)C=C1)C1=CC=C(F)C=C1>>CC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[N+]([O-])=O |f:1.2|"
    )
    response2 = client.post("/autoupdate_reaction_table", json=data2)

    assert response2.status_code == 200
    assert (
        response2.json["novelCompound"] is True
        and "Reactant 2 not in database" in response2.json["reactionTable"]
    )

    # replace first compound with a db compound
    data3 = get_test_data(
        reaction_smiles="C.CC>>CC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[N+]([O-])=O |f:1.2|"
    )
    response3 = client.post("/autoupdate_reaction_table", json=data3)

    assert response3.status_code == 200
    assert (
        response3.json["novelCompound"] is True
        and "Product 1 not in database" in response3.json["reactionTable"]
    )


def test_demo(client: FlaskClient):
    """Tests the reaction table /_process response when not logged in"""
    # no need to login in demo mode
    data = get_test_data(reaction_smiles="C.CC>>CP.CCP", demo="demo")

    response = client.get("/autoupdate_reaction_table", json=data)

    assert_reaction_table_response_for_test_compounds(response)

    # check it still works after login
    login(client)
    response = client.get("/autoupdate_reaction_table", json=data)
    assert_reaction_table_response_for_test_compounds(response)


# polymer fixes coming in later PR
# def test_autoupdate_reaction_table_compounds_polymer_reactions(client: FlaskClient):
#     login(client)
#
#     data = get_test_data(reaction_smiles="PC{-}(CCC{+n}(C)P)C>>C |$;;;;;;;;*;;;*;$|")
#
#     response = client.post("/autoupdate_reaction_table", json=data)
#
#     assert response.status_code == 200
#


def get_test_data(
    reaction_smiles: str,
    demo: str = "not%20demo",
    tutorial: str = "no",
    workgroup: str = "Test-Workgroup",
    workbook: str = "Test-Workbook",
    reaction_id: str = "TW1-001",
    polymer: str = "false",
    polymer_indices=None,
):
    """
    Mimics front end json data that will be sent to autosubmit route
    """
    if polymer_indices is None:
        polymer_indices = []
    return {
        "reaction_smiles": reaction_smiles,
        "demo": demo,
        "tutorial": tutorial,
        "workgroup": workgroup,
        "workbook": workbook,
        "reaction_id": reaction_id,
        "polymer": polymer,
        "polymer_indices": polymer_indices,
    }


def make_url(
    reactants: str,
    products: str,
    demo: str = "not%20demo",
    workgroup: str = "Test-Workgroup",
    reaction_smiles: Optional[str] = None,
    workbook: str = "Test-Workbook",
    reaction_id: str = "TW1-001",
    polymer: str = "false",
    polymerIndices: str = "",
) -> str:
    """
    Constructs a URL for making a GET request based on reactants, products, demo, workgroup, workbook, and reaction ID.

    Args:
        reactants (str): The reactants of the chemical reaction, separated by commas.
        products (str): The products of the chemical reaction, separated by commas.
        demo (str, optional): Whether it's a demo or not. Defaults to "not%20demo".
        workgroup (str, optional): The workgroup associated with the reaction. Defaults to "Test-Workgroup".
        reaction_smiles (str, optional): The SMILES string for the reaction. Defaults to None.
        workbook (str, optional): The workbook associated with the reaction. Defaults to "Test-Workbook".
        reaction_id (str, optional): The ID of the reaction. Defaults to "TW1-001".
        polymer (str, optional): Whether polymer mode is on. Defaults to "false".
        polymerIndices (str, optional): Indices of polymers present. Defaults to "".

    Returns:
        str: The constructed URL for making a GET request.
    """
    # Constructing reaction_smiles from reactants and products
    reactants_list = reactants.split(",")
    products_list = products.split(",")
    if reaction_smiles is None:
        reaction_smiles = ".".join(reactants_list) + ">>" + ".".join(products_list)

    # Defining the endpoint
    endpoint = "/_process"

    # Constructing the URL with the arguments
    url = f"{endpoint}?reactants={reactants}&products={products}&reactionSmiles={reaction_smiles}&demo={demo}&workgroup={workgroup}&workbook={workbook}&reaction_id={reaction_id}&polymer={polymer}&polymerIndices={polymerIndices}"

    return url


def assert_reaction_table_response_for_test_compounds(response):
    assert response.status_code == 200

    soup = BeautifulSoup(response.json["reactionTable"], "html.parser")

    # test names and mol wts
    assert soup.find("input", id="js-reactant1")["value"] == "Testoic Acid"
    assert soup.find("input", id="js-reactant-molecular-weight1")["value"] == "123.0"

    assert soup.find("input", id="js-reactant2")["value"] == "Testamine"
    assert soup.find("input", id="js-reactant-molecular-weight2")["value"] == "101.0"

    assert soup.find("input", id="js-product1")["value"] == "Testproductine"
    assert soup.find("input", id="js-product-molecular-weight1")["value"] == "225.0"

    assert soup.find("input", id="js-product2")["value"] == "Testproductone"
    assert soup.find("input", id="js-product-molecular-weight2")["value"] == "255.0"

    # test hazard codes
    assert (
        soup.find("textarea", id="js-reactant-hazards1").get_text(strip=True)
        == "H900-H901"
    )
    assert (
        soup.find("textarea", id="js-reactant-hazards2").get_text(strip=True) == "H900"
    )

    assert soup.find("textarea", id="js-product-hazard1").get_text(strip=True) == "H900"
    assert soup.find("textarea", id="js-product-hazard2").get_text(strip=True) == "H900"

    # assert table numbers
    assert soup.find("input", id="js-reactant-table-number1")["value"] == "1"
    assert soup.find("input", id="js-reactant-table-number2")["value"] == "2"

    assert soup.find("input", id="js-product-table-number1")["value"] == "3"
    assert soup.find("input", id="js-product-table-number2")["value"] == "4"
