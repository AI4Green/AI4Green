from flask.testing import FlaskClient
from tests.utils import login


def test_autosave(client: FlaskClient):
    """Tests the response when autosaving a reaction"""
    login(client)
    # Define form data to send in the request
    form_data = make_autosave_form()
    # Send a POST request to the route with the form data and confirm response
    response = client.post("/_autosave", data=form_data)
    assert (
        response.status_code == 200 and response.json["feedback"] == "Reaction Updated!"
    )
    assert "need to test time of update" == float


def make_autosave_form():
    return {
        "reactionName": "test reaction name",
        "reactionID": "TW1-001",
        "reactionDescription": "testing routes and services",
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "reactionSmiles": "C.CC>>CP.CCP",
        "userEmail": "{{ current_user.email }}",
        "reactantPrimaryKeys": "1;2;",
        "productPrimaryKeys": "6;7",
        "reagentPrimaryKeys": "3;4",
        "solventPrimaryKeys": "5",
        "numberOfReactants": "2",
        "limitingReactantTableNumber": "1",
        "reactantMasses": "1;0.37;",
        "reactantMassesRaw": "1;0.3691451031772027;",
        "reactantAmounts": "0.01;0.01;",
        "reactantAmountsRaw": "0.008188666885031117;0.008188666885031117;",
        "reactantVolumes": "0.00;0.00;",
        "reactantVolumesRaw": "0.0005285714285714286;",
        "reactantEquivalents": "1;1;",
        "reactantPhysicalForms": "1;2;",
        "reactantDensities": "1.316;0.7;",
        "reactantConcentrations": ";;",
        "reagentNames": "",
        "reagentDensities": "",
        "reagentConcentrations": "",
        "reagentMolecularWeights": "",
        "reagentAmounts": "",
        "reagentAmountsRaw": "",
        "reagentEquivalents": "",
        "reagentPhysicalForms": "",
        "reagentHazards": "",
        "reagentMasses": "",
        "reagentMassesRaw": "",
        "reagentVolumes": "",
        "reagentVolumesRaw": "",
        "reagentSmiles": "",
        "solventPhysicalForms": "",
        "solventNames": "",
        "solventConcentrations": "",
        "solventVolumes": "",
        "solventHazards": "",
        "productPhysicalForms": "1;3;",
        "productAmounts": "0.01;0.01;",
        "productAmountsRaw": "0.008188666885031117;0.008188666885031117;",
        "productMasses": "1.22;0.15;",
        "productMassesRaw": "1.2216672125777923;0.14751883393383558;",
        "mainProductTableNumber": "1",
        "amountUnits": "mmol",
        "massUnits": "mg",
        "volumeUnits": "mL",
        "solventVolumeUnits": "mL",
        "productAmountUnits": "mmol",
        "productMassUnits": "mg",
        "realProductMass": "",
        "unreactedReactantMass": "",
        "reactionTemperature": "",
        "elementSustainability": "undefined",
        "batchFlow": "-select-",
        "isolationMethod": "undefined",
        "catalystUsed": "-select-",
        "catalystRecovered": "-select-",
        "otherHazardTextArea": "",
        "customProtocol1": "",
        "customProtocol2": "",
        "selectedRadioButtons": "",
        "researcher": "",
        "supervisor": "",
        "complete": "not complete",
        "reactantNames": "Benzoic Acid;Ethylamine;",
        "reactantMolecularWeights": "122.12;45.08;",
        "reactantHazards": "H315-H318-H372;H220-H319-H335;",
        "reactantPhysicalFormsText": "Dense solid;Non-volatile liquid;",
        "reagentPhysicalFormsText": "",
        "solventPhysicalFormsText": "",
        "productNames": "N-Ethylbenzamide;Water;",
        "productMolecularWeights": "149.19;18.015;",
        "productHazards": "H302;Unknown;",
        "productPhysicalFormsText": "Dense solid;Unknown;",
        "summary_to_print": "no summary data",
        "massEfficiency": "",
        "conversion": "",
        "selectivity": "",
        "toExport": '[{"key":"Temperature Sustainability"},{"key":"Elements Sustainability"},{"key":"Batch or Flow Sustainability"},{"key":"Isolation Sustainability"},{"key":"Catalyst Sustainability"},{"key":"Recovery Sustainability"},{"key":"Atom Economy Sustainability"},{"key":"Mass Efficiency Sustainability"},{"key":"Yield Sustainability"},{"key":"Conversion Sustainability"},{"key":"Selectivity Sustainability"}]',
    }
