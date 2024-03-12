from flask import Flask
from flask.testing import FlaskClient
from sources import services
from tests.utils import assert_expected_values, login


def test_summary(app: Flask, client: FlaskClient):
    """Tests the response when adding a novel compound from the sketcher"""
    login(client)
    # Define form data to send in the request
    form_data = make_summary_form()
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_summary", data=form_data)
    assert response.status_code == 200 and response.json["ae"] == 16


# def test_summary_demo():
#     pass


def make_summary_form():
    # Provided data dictionary
    return {
        # reaction
        "reactionSmiles": "C.CC>>CP.CCP",
        "reactionName": "test reaction name",
        "reactionID": "TW1-001",
        "reactionDescription": "testing routes and services",
        "demo": "not demo",
        "tutorial": "no",
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        # reactant
        "reactants": "Testoic Acid;Testamine",
        "reactantIds": "-1;-2",
        "reactantMolecularWeights": "123;101",
        "roundedReactantMasses": "123;101",
        "reactantMasses": "123;101",
        "roundedReactantAmounts": "1.00;10.0",
        "reactantAmounts": "1.00;10.0",
        "roundedReactantVolumes": "0.10;0.02",
        "reactantDensities": "1.204;0.995",
        "reactantVolumes": "0.10;0.02",
        "reactantEquivalents": "1;10",
        "reactantConcentrations": "-;-",
        "reactantPhysicalForms": "Dense solid;Dense solid",
        "limitingReactantTableNumber": "1",
        "reactantHazards": "H900-H901;H900",
        "reactantMassSum": "360",
        "reactantMolecularWeightSum": "270",
        "reactantPrimaryKeys": "1;2",
        # reagents
        "reagentIds": "-3;-4",
        "reagentNames": "Testide;Testol",
        "roundedReagentAmounts": "0.20;2.60",
        "reagentAmounts": "0.20;2.60",
        "reagentEquivalents": "0.2;2.6",
        "reagentMolecularWeights": "150;175",
        "reagentDensities": ";0.703",
        "reagentConcentrations": ";;",
        "roundedReagentMasses": "3.21;297",
        "reagentMasses": "3.2054;296.706",
        "roundedReagentVolumes": ";0.42",
        "reagentVolumes": ";0.422475",
        "reagentPhysicalForms": "Dense solid;Highly volatile liquid",
        "reagentHazards": "H901;H901",
        "reagentMolecularWeightSum": "361",
        "reagentPrimaryKeys": "3;4",
        "reagentTableNumbers": "3;4",
        "reagents": "Testide;Testol",
        # solvents
        # 'solventNames': ['Testanol'],
        "solventConcentrations": "0.50",
        "solventIds": "-5",
        "solventVolumes": "2",
        "solventPhysicalForms": "Volatile liquid",
        "solventHazards": "H900",
        "solventTableNumbers": "5",
        "solvents": "Testanol",
        "solventPrimaryKeys": "5",
        "numberOfSolvents": "1",
        # products
        "productIds": "-6;-7",
        "products": "Testproductine;Testproductone",
        "roundedProductMasses": "93.0;18.0",
        "productMasses": "93.0;18.0",
        "roundedProductAmounts": "1.00;1.00",
        "productAmounts": "1.00;1.00",
        "productMolecularWeights": "225;255",
        "productHazards": "H900;Unknown",
        "productPrimaryKeys": "6;7",
        "productPhysicalForms": "Unknown;Unknown",
        "mainProductTableNumber": "7",
        "productTableNumbers": "6;7",
        "mainProduct": "1",
        # units
        "productMassUnit": "mg",
        "massUnit": "mg",
        "solventVolumeUnit": "mL",
        "volumeUnit": "mL",
        "amountUnit": "mmol",
        # summary data
        "js_summary_table_data": "no data",
        "realProductMass": "90",
        "unreactedReactantMass": "2",
        "reactionTemperature": "70",
        "elementSustainability": "3",
        "batchFlow": "Batch",
    }
