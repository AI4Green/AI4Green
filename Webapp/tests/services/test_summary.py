from flask.testing import FlaskClient
from tests.utils import login


def test_summary(client: FlaskClient):
    """Tests the response when loading the summary table for a workbook reaction"""
    login(client)
    # Define form data to send in the request
    form_data = make_workbook_summary_form()
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_summary", data=form_data)
    assert response.status_code == 200
    assert_summary_table_html_formed(response)


def test_summary_demo(client: FlaskClient):
    """Tests the response when loading the summary table for a demo reaction"""
    # Define form data to send in the request
    form_data = make_demo_summary_form()
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_summary", data=form_data)
    assert response.status_code == 200
    assert_summary_table_html_formed(response)


def test_summary_tutorial(client: FlaskClient):
    """Tests the response when loading the summary table for a tutorial reaction"""
    # Define form data to send in the request
    form_data = make_tutorial_summary_form()
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_summary", data=form_data)
    assert response.status_code == 200
    assert_summary_table_html_formed(response)


def assert_summary_table_html_formed(response):
    # element sustainability has been calculated
    assert (
        """<td id="js-elements-cell" class="hazard-acceptable"><select size="1" id="js-elements" class="hazard-acceptable">"""
        in response.json["summary"]
    ), "Element sustainability has not been calculated"
    # atom economy
    assert (
        """<td id="js-ae-cell" class="hazard-hazardous"><input id="js-ae" type="number" class="hazard-hazardous to-export" name="Atom Efficiency" value="40.4" style="width: 80px; border:none;" readonly></td>"""
        in response.json["summary"]
    ), "Atom economy has not been calculated"
    # theoretical yield
    assert (
        """<td colspan="2"><input type="number" id="js-product-rounded-mass1"
                               value="93.0" style="width:80px; border:none" readonly></td>"""
        in response.json["summary"]
    ), "Theoretical yield has not been calculated"


def make_demo_summary_form():
    summary_form = make_default_summary_form()
    summary_form.update(
        {
            "demo": "demo",
            "tutorial": "no",
        }
    )
    return summary_form


def make_tutorial_summary_form():
    summary_form = make_default_summary_form()
    summary_form.update(
        {
            "demo": "not demo",
            "tutorial": "yes",
        }
    )
    return summary_form


def make_workbook_summary_form():
    """All the fields that are sent from JS to /_summary endpoint"""
    summary_form = make_default_summary_form()
    summary_form.update(
        {
            "reactionName": "test reaction name",
            "reactionID": "TW1-001",
            "reactionDescription": "testing routes and services",
            "demo": "not demo",
            "tutorial": "no",
            "workgroup": "Test-Workgroup",
            "workbook": "Test-Workbook",
        }
    )
    return summary_form


def make_default_summary_form():
    return {
        # reaction
        "reactionSmiles": "C.CC>>CP.CCP",
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
        "productHazards": "H900;H900",
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
