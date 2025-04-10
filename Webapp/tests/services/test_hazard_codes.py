from sources import services


def test_toxicity_warnings():
    """Tests the toxicity warning functions return the correct toxicity type for the hazard codes"""
    hazard_code_list = []
    assert services.hazard_code.get_toxicities(hazard_code_list) == []
    hazard_code_list = ["H220", "H224"]
    assert services.hazard_code.get_toxicities(hazard_code_list) == []
    hazard_code_list = ["H220", "H350"]
    assert services.hazard_code.get_toxicities(hazard_code_list) == ["carcinogen"]
    hazard_code_list = ["H340", "H220"]
    assert services.hazard_code.get_toxicities(hazard_code_list) == ["mutagen"]
    hazard_code_list = ["H360", "H221", "H226", "H220"]
    assert services.hazard_code.get_toxicities(hazard_code_list) == [
        "reproductive toxin"
    ]
    hazard_code_list = ["H360", "H221", "H350", "H220", "H340"]
    assert services.hazard_code.get_toxicities(hazard_code_list) == [
        "carcinogen",
        "mutagen",
        "reproductive toxin",
    ]
