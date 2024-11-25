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


def make_autosave_form():
    return {
        "reactionName": "test reaction name",
        "reactionID": "TW1-001",
        "reactionDescription": "testing routes and services",
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "reactionSmiles": "C.CC>>CP.CCP",
        "reactionRXN": exampleRXN(),
        "reactionImage": "data:image/jpeg;base64,/9j/4AAQSkZJRgABAQAAAQABAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCADrAlgDAREAAhEBAxEB/8QAHwAAAQUBAQEBAQEAAAAAAAAAAAECAwQFBgcICQoL/8QAtRAAAgEDAwIEAwUFBAQAAAF9AQIDAAQRBRIhMUEGE1FhByJxFDKBkaEII0KxwRVS0fAkM2JyggkKFhcYGRolJicoKSo0NTY3ODk6Q0RFRkdISUpTVFVWV1hZWmNkZWZnaGlqc3R1dnd4eXqDhIWGh4iJipKTlJWWl5iZmqKjpKWmp6ipqrKztLW2t7i5usLDxMXGx8jJytLT1NXW19jZ2uHi4+Tl5ufo6erx8vP09fb3+Pn6/8QAHwEAAwEBAQEBAQEBAQAAAAAAAAECAwQFBgcICQoL/8QAtREAAgECBAQDBAcFBAQAAQJ3AAECAxEEBSExBhJBUQdhcRMiMoEIFEKRobHBCSMzUvAVYnLRChYkNOEl8RcYGRomJygpKjU2Nzg5OkNERUZHSElKU1RVVldYWVpjZGVmZ2hpanN0dXZ3eHl6goOEhYaHiImKkpOUlZaXmJmaoqOkpaanqKmqsrO0tba3uLm6wsPExcbHyMnK0tPU1dbX2Nna4uPk5ebn6Onq8vP09fb3+Pn6/9oADAMBAAIRAxEAPwD+/igAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgDBj8VeGJfE1z4Ki8R6DJ4ys9EtPE154Sj1jT38TWvhu/vrvTLHxBc6CtwdVg0S91KwvtPtNVltFsLi+s7u0huHnt5o0AIz4v8JDxaPAJ8UeHR46bw63i9fBR1vTP+Etbwkmpror+KB4c+0/2wfDq6w6aS2tiz/sxdTdbA3IumERAOioA8O8c/tO/s1/DDx94b+FPxK/aG+Bvw8+KPjJbZvCHw28c/FnwF4S8feKlvLh7Szbw34O1/X9P8Q64t3dRvbWx0vTroT3CPDFvkUqAD3EEEAgggjII5BB6EHuDQAUAFABQAUAFABQAUABIAJJAAGSTwAB1JPYCjbcN9jkdB+IHgPxTdPY+GPG3hHxHexxGaSz0HxJo2r3UcIIUyvb6fe3EqxBmVTIyBQWAJyRXiZdxLw5m9aWHynP8kzTERg6kqGXZrgcbWjBOznKlhq9Wagm0nJxsm7XPXx/D+fZXSVfM8kzfLqMpKEa2Py3G4OlKbu1FVMRRpwcmk2op3dnpodBHqemy6jcaRFqFjJq1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fDejaZq/jzXYNR8T+Jb+HV/FN5quoRXviHXLpLhZ9Wv5LgA9coAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgDx74b/s8/AD4N+Kvib46+EPwN+D3wq8b/GvX18V/GXxj8N/hn4L8D+Kvi14oXUdf1hfEnxN8Q+GNE0vV/Hmvrq/irxRqi6x4qvNV1Eaj4k1++FyLnWNRluQD2GgAoAKACgAoAKACgAoAKACgD82f2rv+CYHwM/aN+IVr+0d8PPEXjr9kX9tbQbSKDwx+2J+zTqNn4K+KWoQ2DaFcab4U+NOjSWlz4J/aN+Fclz4U8LWmu/DT4y6D4n0vUPDekN4a0i/8OWd/dzOAfLY1b/gvvcN/wAMuP4V/YssNQW0d5v+CqMOsatdeFJPBhujpEIsf2DpSviZP2nWsiPE89nefEV/2bLa7Rlj1i8gkTwsoB9W/sm/8ExvgL+zP491L9oXxfq/jb9qn9s/xTZi28Z/tj/tJ6pb+PfjE0Mw103fhr4WQyW0HhP9n74ZxjxPrumaZ8NvgxofhDQY/Dc+naDrj+IodG065hAP0doAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoAKACgAoA//Z",
        "polymerMode": "",
        "polymerIndices": "[]",
        "polymerisationType": "",
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
        "polymerMn": "",
        "polymerMw": "",
        "polymerDispersity": "",
        "polymerMassMethod": "",
        "polymerMassCalibration": "",
        "polymerTg": "",
        "polymerTm": "",
        "polymerTc": "",
        "polymerThermalMethod": "",
        "polymerThermalCalibration": "",
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
        "productIntendedDPs": "",
        "summary_to_print": "no summary data",
        "massEfficiency": "",
        "conversion": "",
        "selectivity": "",
        "toExport": '[{"key":"Temperature Sustainability"},{"key":"Elements Sustainability"},{"key":"Batch or Flow Sustainability"},{"key":"Isolation Sustainability"},{"key":"Catalyst Sustainability"},{"key":"Recovery Sustainability"},{"key":"Atom Economy Sustainability"},{"key":"Mass Efficiency Sustainability"},{"key":"Yield Sustainability"},{"key":"Conversion Sustainability"},{"key":"Selectivity Sustainability"}]',
    }


def exampleRXN():
    return """$RXN

             -INDIGO- 1111241325

              2  2
            $MOL

              -INDIGO-11112413252D

              1  0  0  0  0  0  0  0  0  0999 V2000
                3.0500   -4.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            M  END
            $MOL

              -INDIGO-11112413252D

              2  1  0  0  0  0  0  0  0  0999 V2000
                6.3500   -4.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                7.3500   -4.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
              1  2  1  0  0  0  0
            M  END
            $MOL

              -INDIGO-11112413252D

              2  1  0  0  0  0  0  0  0  0999 V2000
               15.9000   -4.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               16.9000   -4.7000    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
              1  2  1  0  0  0  0
            M  END
            $MOL

              -INDIGO-11112413252D

              3  2  0  0  0  0  0  0  0  0999 V2000
               19.4000   -4.6750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               20.4000   -4.6750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               20.9000   -5.5410    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
              1  2  1  0  0  0  0
              2  3  1  0  0  0  0
            M  END"""
