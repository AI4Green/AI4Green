from flask.testing import FlaskClient
from tests.utils import assert_expected_values, login


def test_reaction_table_compounds_present(client: FlaskClient):
    login(client)
    response = client.get(
        "/_process?reactants=C,CC&products=CP,CCP&reactionSmiles=C.CC>>CP.CCP&demo=not%20demo&workgroup=Test-Workgroup&workbook=Test-Workbook&reaction_id=TW1-001"
    )
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
        """<td><textarea class="text-in-table" id="js-reactant-hazards2" readonly>H789</textarea></td>"""
        in response.json["reactionTable"]
    )
    assert (
        """<td><input type="number" id="js-product-table-number2" value="4" style="width:40px; border:none" readonly></td>"""
        in response.json["reactionTable"]
    )
    # assert response.json['novelCompound'] is True and 'Reactant 1 not in database' in response.json['reactionTable']


def test_reaction_table_ion_reaction(client: FlaskClient):
    login(client)
    # response = client.get()


#
#
# def test_reaction_table_compounds_not_present(client: FlaskClient):
#     """Tests the response when using reactants and products not in the database outside demo/tutorial mode"""
#     login(client)
#     # Define GET URL # todo - break down into pieces
#
#     # Send a POST request to the route with the form data
#     # Although this is the default example reaction, these are not present in the testing database (yet at least)
#     response = client.get("/_process?reactants=OC(=O)C1=CC=CC=C1,CCN&products=CCNC(=O)C1=CC=CC=C1&reactionSmiles=OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1&demo=not%20demo&workgroup=Test-Workgroup&workbook=Test-Workbook&reaction_id=TW1-001")
#     assert response.status_code == 200
#     assert response.json['novelCompound'] is True and 'Reactant 1 not in database' in response.json['reactionTable']
#
