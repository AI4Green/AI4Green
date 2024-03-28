from typing import Optional

from flask.testing import FlaskClient
from tests.utils import assert_expected_values, login


def test_reaction_table_compounds_present(client: FlaskClient):
    """
    Tests the reaction table /_process when the compounds used are present in the Compound table of the database
    We test the results by looking for substrings of the HTML which include the compound data
    """
    login(client)

    url = make_url(reactants="C,CC", products="CP,CCP")

    response = client.get(url)
    assert_reaction_table_response_for_test_compounds(response)


def test_reaction_table_ion_reactions(client: FlaskClient):
    """
    Tests the reaction table /_process when the compounds used are ions not present in the database
    """
    login(client)
    # with all ions in
    url = make_url(
        reactants="CC(C(=O)N1CCCC1C(=O)OC(C)(C)C)[Nplus]([Ominus])=O,F[Bminus](F)(F)F,FC1=CC=C(C=C1)[Asplus](C1=CC=C(F)C=C1)(C1=CC=C(F)C=C1)C1=CC=C(F)C=C1",
        products="CC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[Nplus]([Ominus])=O",
        reaction_smiles="CC(C(=O)N1CCCC1C(=O)OC(C)(C)C)[Nplus]([Ominus])=O.F[Bminus](F)(F)F.FC1=CC=C(C=C1)[Biplus](C1=CC=C(F)C=C1)(C1=CC=C(F)C=C1)C1=CC=C(F)C=C1%3E%3ECC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[Nplus]([Ominus])=O%20|f:1.2|",
    )
    response = client.get(url)
    assert response.status_code == 200
    assert (
        response.json["novelCompound"] is True
        and "Reactant 1 not in database" in response.json["reactionTable"]
    )

    # replace first ion with database compound
    url = make_url(
        reactants="C,F[Bminus](F)(F)F,FC1=CC=C(C=C1)[Asplus](C1=CC=C(F)C=C1)(C1=CC=C(F)C=C1)C1=CC=C(F)C=C1",
        products="CC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[Nplus]([Ominus])=O",
        reaction_smiles="C.F[Bminus](F)(F)F.FC1=CC=C(C=C1)[Biplus](C1=CC=C(F)C=C1)(C1=CC=C(F)C=C1)C1=CC=C(F)C=C1%3E%3ECC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[Nplus]([Ominus])=O%20|f:1.2|",
    )
    response = client.get(url)
    assert response.status_code == 200
    assert (
        response.json["novelCompound"] is True
        and "Reactant 2 not in database" in response.json["reactionTable"]
    )

    # with chiral product
    url = make_url(
        reactants="C",
        products="CC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[Nplus]([Ominus])=O",
        reaction_smiles="C%3E%3ECC(C)(C)OC(=O)C1CCCN1C(=O)[C@@](C)(C1=CC=C(F)C=C1)[Nplus]([Ominus])=O",
    )
    response = client.get(url)
    assert response.status_code == 200
    assert (
        response.json["novelCompound"] is True
        and "Product 1 not in database" in response.json["reactionTable"]
    )


def test_reaction_table_compounds_not_present(client: FlaskClient):
    """Tests the reaction table /_process response when using neutral reactants and products not in the database"""
    login(client)
    # Although this is the default example reaction, these are not in db. Defined in - (tests/populate_database.py)
    url = make_url(reactants="OC(=O)C1=CC=CC=C1,CCN", products="CCNC(=O)C1=CC=CC=C1")
    response = client.get(url)
    assert response.status_code == 200
    assert (
        response.json["novelCompound"] is True
        and "Reactant 1 not in database" in response.json["reactionTable"]
    )


def test_demo(client: FlaskClient):
    """Tests the reaction table /_process response when not logged in"""
    # no need to login in demo mode
    url = make_url(
        reactants="C,CC",
        products="CP,CCP",
        demo="demo",
    )
    response = client.get(url)
    assert_reaction_table_response_for_test_compounds(response)
    # check it still works after login
    login(client)
    url = make_url(
        reactants="C,CC",
        products="CP,CCP",
        demo="demo",
    )
    response = client.get(url)
    assert_reaction_table_response_for_test_compounds(response)


def make_url(
    reactants: str,
    products: str,
    demo: str = "not%20demo",
    workgroup: str = "Test-Workgroup",
    reaction_smiles: Optional[str] = None,
    workbook: str = "Test-Workbook",
    reaction_id: str = "TW1-001",
) -> str:
    """
    Constructs a URL for making a GET request based on reactants, products, demo, workgroup, workbook, and reaction ID.

    Args:
        reactants (str): The reactants of the chemical reaction, separated by commas.
        products (str): The products of the chemical reaction, separated by commas.
        demo (str, optional): Whether it's a demo or not. Defaults to "not%20demo".
        workgroup (str, optional): The workgroup associated with the reaction. Defaults to "Test-Workgroup".
        workbook (str, optional): The workbook associated with the reaction. Defaults to "Test-Workbook".
        reaction_id (str, optional): The ID of the reaction. Defaults to "TW1-001".

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
    url = f"{endpoint}?reactants={reactants}&products={products}&reactionSmiles={reaction_smiles}&demo={demo}&workgroup={workgroup}&workbook={workbook}&reaction_id={reaction_id}"

    return url


def assert_reaction_table_response_for_test_compounds(response):
    assert response.status_code == 200
    assert (
        """<td><input style="width: 70px; border-width:0px;" type="number" step="any" value="123.0"
                    id="js-reactant-molecular-weight1" readonly>
            </td>"""
        in response.json["reactionTable"]
    )
    assert (
        """<td style='width: 200px;'><input  id="js-product1" style="width: 200px; border-width:0px;" value="Testproductine"></td>"""
        in response.json["reactionTable"]
    )
    assert (
        """<td><textarea class="text-in-table" id="js-reactant-hazards2" readonly>H900</textarea></td>"""
        in response.json["reactionTable"]
    )
    assert (
        """<td><input type="number" id="js-product-table-number2" value="4" style="width:40px; border:none" readonly></td>"""
        in response.json["reactionTable"]
    )
