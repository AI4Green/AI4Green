"""
Tests for reaction workups and purifications (/api/workups).

Sets up its own workgroup so it doesn't affect other test modules:
- Workup-Workgroup / Workup-Workbook: the main test user ("test_username") and
  "wu_member" are members. The main test user created WU1-001 and WU1-003
  (editable) and WU1-002 (locked).
- Workup-Other: a second workbook in the same workgroup, for cross-workbook checks.
- Workup-Private: a workbook only "wu_member" is in, holding a novel compound the
  main test user mustn't be able to link.
- "wu_outsider" is in no workbook here, to check non-members get 404s.

Materials use the seeded compound database: "Testanol" is a Compound linked to the
CHEM21 solvent "Testanol" (flag 2). The fixture adds "Workupane", a CHEM21 solvent
with no Compound row, and novel compounds "Workupol" (Workup-Workbook) and
"Secretol" (Workup-Private).
"""

from datetime import datetime
from decimal import Decimal
from typing import Any, Dict, Iterator, List

import pytest
from flask import Flask
from flask.testing import FlaskClient
from sources import models, services
from sources.extensions import db
from sqlalchemy.exc import IntegrityError
from tests.populate_database import add_test_user
from tests.utils import add_workgroup, login, main_test_person

WORKGROUP = "Workup-Workgroup"
WORKBOOK = "Workup-Workbook"
OTHER_WORKBOOK = "Workup-Other"
PASSWORD = "wu_pw"
TEMPLATES = f"/api/workups/{WORKGROUP}/{WORKBOOK}/templates"


def _add_reaction(
    person: models.Person, workbook: models.WorkBook, code: str, complete: str
) -> models.Reaction:
    """Adds a minimal reaction to the workbook."""
    return models.Reaction.create(
        creator=person.id,
        time_of_creation=datetime.now(),
        reaction_id=code,
        name=f"workup reaction {code}",
        workbooks=workbook.id,
        reaction_table_data={},
        summary_table_data={},
        complete=complete,
        status="active",
    )


@pytest.fixture(scope="module")
def workup_setup(app: Flask) -> Dict[str, int]:
    """
    Creates the users, workbooks and reactions for these tests.

    Returns:
        Ids the tests need.
    """
    with app.app_context():
        owner = main_test_person()
        member = add_test_user(
            "wu_member", "wu_member@test.com", "Workup Member", PASSWORD, True
        )
        add_test_user(
            "wu_outsider", "wu_outsider@test.com", "Workup Outsider", PASSWORD, True
        )

        workgroup = add_workgroup(
            WORKGROUP, principal_investigator=[owner], standard_member=[member]
        )
        # no shared helper for workbooks yet; mirrors tests/populate_database.py
        workbook = models.WorkBook.create(
            name=WORKBOOK, group=workgroup.id, users=[owner, member], abbreviation="WU1"
        )
        other_workbook = models.WorkBook.create(
            name=OTHER_WORKBOOK, group=workgroup.id, users=[owner], abbreviation="WU2"
        )
        private_workbook = models.WorkBook.create(
            name="Workup-Private",
            group=workgroup.id,
            users=[member],
            abbreviation="WU3",
        )
        models.Solvent.create(name="Workupane", flag=4, hazard="H900")
        workupol = models.NovelCompound.create(
            name="Workupol",
            smiles="CCCCCCCO",
            molec_formula="C7H16O",
            workbook=workbook.id,
        )
        secretol = models.NovelCompound.create(
            name="Secretol",
            smiles="CCCCCCCCO",
            molec_formula="C8H18O",
            workbook=private_workbook.id,
        )
        _add_reaction(owner, workbook, "WU1-001", "not complete")
        _add_reaction(owner, workbook, "WU1-002", "complete")
        _add_reaction(owner, workbook, "WU1-003", "not complete")
        db.session.commit()
        return {
            "workbook_id": workbook.id,
            "other_workbook_id": other_workbook.id,
            "testanol_id": services.compound.from_name("Testanol").id,
            "workupol_id": workupol.id,
            "secretol_id": secretol.id,
        }


@pytest.fixture(autouse=True)
def clean_workups(app: Flask, workup_setup: Dict[str, int]) -> Iterator[None]:
    """Removes every workup step and template after each test."""
    yield
    with app.app_context():
        db.session.query(models.WorkupStepMaterial).delete()
        db.session.query(models.WorkupStep).delete()
        db.session.query(models.WorkupTemplateStep).delete()
        db.session.query(models.WorkupTemplate).delete()
        db.session.commit()


@pytest.fixture(autouse=True)
def audit_messages(monkeypatch: pytest.MonkeyPatch) -> List[Any]:
    """Captures workup audit messages instead of sending them."""
    sent: List[Any] = []
    monkeypatch.setattr(services.reaction_editing_history, "send_message", sent.append)
    return sent


def _url(code: str = "WU1-001", workbook: str = WORKBOOK) -> str:
    """The workup URL for one reaction."""
    return f"/api/workups/{WORKGROUP}/{workbook}/{code}"


def _add_step(
    client: FlaskClient,
    step_type: str,
    category: str = "WORKUP",
    code: str = "WU1-001",
) -> Dict[str, Any]:
    """Adds a step and returns it."""
    response = client.post(
        f"{_url(code)}/steps", json={"category": category, "stepType": step_type}
    )
    assert response.status_code == 201, response.json
    return response.json


def _steps(
    client: FlaskClient, code: str = "WU1-001", category: str = "WORKUP"
) -> Dict[str, Any]:
    """Gets a reaction's steps for one category."""
    response = client.get(_url(code), query_string={"category": category})
    assert response.status_code == 200, response.json
    return response.json


def _save(
    client: FlaskClient, step: Dict[str, Any], code: str = "WU1-001", **body: Any
) -> Any:
    """Saves a step with the given body fields (data defaults to {})."""
    return client.put(f"{_url(code)}/steps/{step['id']}", json={"data": {}, **body})


def _field_names(messages: List[Any]) -> List[str]:
    """The field_name of each captured audit message."""
    return [message.field_name for message in messages]


# --- Steps --------------------------------------------------------------------


def test_creator_adds_and_lists_steps(
    client: FlaskClient, audit_messages: List[Any]
) -> None:
    """Steps are appended in order and only listed under their own category."""
    login(client)
    wash = _add_step(client, "WASH")
    extraction = _add_step(client, "EXTRACTION")
    _add_step(client, "FLASH_CHROMATOGRAPHY", category="PURIFICATION")

    assert (wash["sortOrder"], extraction["sortOrder"]) == (1, 2)
    workup = _steps(client)
    assert workup["canEdit"] is True
    assert [step["stepType"] for step in workup["steps"]] == ["WASH", "EXTRACTION"]
    purification = _steps(client, category="PURIFICATION")
    assert [step["stepType"] for step in purification["steps"]] == [
        "FLASH_CHROMATOGRAPHY"
    ]
    assert _field_names(audit_messages) == ["Workup Step Added"] * 3
    assert audit_messages[0].change_details == {
        "category": "WORKUP",
        "step_type": "WASH",
    }


def test_unknown_category_or_step_type_is_rejected(client: FlaskClient) -> None:
    """Category and step type must be known values."""
    login(client)
    bad_type = client.post(
        f"{_url()}/steps", json={"category": "WORKUP", "stepType": "TELEPORT"}
    )
    bad_category = client.post(
        f"{_url()}/steps", json={"category": "LUNCH", "stepType": "WASH"}
    )
    missing_category = client.get(_url())

    assert bad_type.status_code == 400
    assert bad_category.status_code == 400
    assert missing_category.status_code == 400


def test_step_values_are_cleaned_and_saved(client: FlaskClient) -> None:
    """Values are trimmed and typed, empty values and unknown keys are dropped."""
    login(client)
    step = _add_step(client, "WASH")

    response = client.put(
        f"{_url()}/steps/{step['id']}",
        json={
            "data": {"repeats": "3", "not_a_field": "ignored"},
            "notes": " Emulsion formed ",
        },
    )

    assert response.status_code == 200, response.json
    assert response.json["data"] == {"repeats": 3}
    assert response.json["notes"] == "Emulsion formed"
    assert _steps(client)["steps"][0]["data"] == response.json["data"]


@pytest.mark.parametrize(
    "step_type, data",
    [
        ("WASH", {"repeats": "-1"}),
        ("WASH", {"repeats": "lots"}),
        ("WASH", {"repeats": True}),
        ("FILTRATION", {"filter_medium": 5}),
        ("FILTRATION", {"filter_medium": "Coffee filter"}),
        ("EXTRACTION", {"phase_kept": "Gas"}),
        ("PH_ADJUST", {"target_ph": 15}),
        ("WASH", ["not", "an", "object"]),
    ],
)
def test_invalid_step_values_are_rejected(
    client: FlaskClient, step_type: str, data: Any
) -> None:
    """Values that don't fit the step type's definition get a 400 with a message."""
    login(client)
    step = _add_step(client, step_type)

    response = client.put(f"{_url()}/steps/{step['id']}", json={"data": data})

    assert response.status_code == 400
    assert response.json["message"]
    assert _steps(client)["steps"][0]["data"] == {}


def test_edits_are_audited_only_when_something_changes(
    client: FlaskClient, audit_messages: List[Any]
) -> None:
    """Autosaving the same values twice records one edit, with the changed keys."""
    login(client)
    step = _add_step(client, "WASH")
    body = {"data": {"repeats": 2}, "notes": ""}

    client.put(f"{_url()}/steps/{step['id']}", json=body)
    client.put(f"{_url()}/steps/{step['id']}", json=body)

    edits = [m for m in audit_messages if m.field_name == "Workup Step Edited"]
    assert len(edits) == 1
    assert edits[0].change_details["changes"] == {
        "data.repeats": {"old_value": None, "new_value": 2}
    }


def test_steps_can_be_reordered_and_deleted(
    client: FlaskClient, audit_messages: List[Any]
) -> None:
    """Reordering needs every step id; deleting closes the gap in the order."""
    login(client)
    ids = [_add_step(client, t)["id"] for t in ("WASH", "EXTRACTION", "FILTRATION")]

    reordered = client.put(
        f"{_url()}/steps/order",
        json={"category": "WORKUP", "stepIds": list(reversed(ids))},
    )
    missing_one = client.put(
        f"{_url()}/steps/order", json={"category": "WORKUP", "stepIds": ids[:2]}
    )
    deleted = client.delete(f"{_url()}/steps/{ids[1]}")

    assert reordered.status_code == 200, reordered.json
    assert [s["stepType"] for s in reordered.json["steps"]] == [
        "FILTRATION",
        "EXTRACTION",
        "WASH",
    ]
    assert missing_one.status_code == 400
    assert deleted.status_code == 200
    remaining = _steps(client)["steps"]
    assert [(s["stepType"], s["sortOrder"]) for s in remaining] == [
        ("FILTRATION", 1),
        ("WASH", 2),
    ]
    assert "Workup Steps Reordered" in _field_names(audit_messages)
    assert "Workup Step Deleted" in _field_names(audit_messages)


def test_step_from_another_reaction_is_not_found(client: FlaskClient) -> None:
    """A step id is only accepted through its own reaction's URL."""
    login(client)
    step = _add_step(client, "WASH", code="WU1-001")

    update = client.put(f"{_url('WU1-003')}/steps/{step['id']}", json={"data": {}})
    delete = client.delete(f"{_url('WU1-003')}/steps/{step['id']}")

    assert update.status_code == 404
    assert delete.status_code == 404


def test_step_definitions_list_each_step_types_fields(client: FlaskClient) -> None:
    """The front end gets the labels, kinds and units for every step type."""
    login(client)
    response = client.get("/api/workups/step-definitions")

    assert response.status_code == 200
    assert set(response.json) == {x.value for x in models.WorkupStepType}
    wash = response.json["WASH"]
    assert [field["key"] for field in wash["fields"]] == ["repeats"]
    assert [slot["slot"] for slot in wash["materials"]] == [
        "organic_phase",
        "aqueous_phase",
    ]
    assert wash["materials"][1]["repeated"] is True
    assert response.json["TEMPERATURE"]["temperature"] == "Target temperature"
    assert response.json["WASH"]["temperature"] is None


# --- Materials and temperature ------------------------------------------------


def test_materials_link_to_compounds_novel_compounds_and_solvents(
    client: FlaskClient, workup_setup: Dict[str, int]
) -> None:
    """
    Each material links to one record and takes its name from it; the CHEM21 flag
    comes from the linked solvent.
    """
    login(client)
    step = _add_step(client, "FLASH_CHROMATOGRAPHY", category="PURIFICATION")

    response = _save(
        client,
        step,
        materials=[
            {
                "slot": "stationary_phase",
                "novelCompoundId": workup_setup["workupol_id"],
                "name": "typed name is replaced",
                "quantity": "40",
                "unit": "g",
            },
            {
                "slot": "eluent",
                "compoundId": workup_setup["testanol_id"],
                "quantity": 300,
            },
            {"slot": "eluent", "solventName": "Workupane", "quantity": 700},
        ],
    )

    assert response.status_code == 200, response.json
    by_name = {m["name"]: m for m in response.json["materials"]}
    assert set(by_name) == {"Workupol", "Testanol", "Workupane"}
    assert by_name["Workupol"]["role"] == "STATIONARY_PHASE"
    assert by_name["Workupol"]["quantity"] == 40
    assert by_name["Workupol"]["flag"] is None
    testanol = by_name["Testanol"]
    assert testanol["role"] == "ELUENT"
    assert testanol["compoundId"] == workup_setup["testanol_id"]
    assert (testanol["unit"], testanol["fate"], testanol["flag"]) == ("mL", "WASTE", 2)
    assert by_name["Workupane"]["solventName"] == "Workupane"
    assert by_name["Workupane"]["flag"] == 4


def test_unlinked_and_empty_materials(client: FlaskClient) -> None:
    """Unknown materials are kept by name; completely empty rows are dropped."""
    login(client)
    step = _add_step(client, "WASH")

    response = _save(
        client,
        step,
        materials=[
            {"slot": "aqueous_phase", "name": " Brine ", "quantity": "20"},
            {"slot": "organic_phase", "name": "", "quantity": ""},
        ],
    )

    assert response.status_code == 200, response.json
    [brine] = response.json["materials"]
    assert brine["name"] == "Brine"
    assert brine["compoundId"] is None and brine["solventName"] is None


@pytest.mark.parametrize(
    "material",
    [
        {"slot": "eluent", "compoundId": 1, "solventName": "Workupane"},
        {"slot": "eluent", "compoundId": 99999999},
        {"slot": "eluent", "compoundId": "1"},
        {"slot": "eluent", "solventName": "Not A Solvent"},
        {"slot": "eluent", "name": "Hexane", "quantity": "-5"},
        {"slot": "eluent", "name": "Hexane", "unit": "g"},
        {"slot": "eluent", "name": "Hexane", "fate": "BURNT"},
        {"slot": "eluent", "name": "x" * 201},
        {"slot": "not_a_slot", "name": "Hexane"},
        "not an object",
    ],
)
def test_invalid_materials_are_rejected(client: FlaskClient, material: Any) -> None:
    """Materials must fit their slot and link to at most one real record."""
    login(client)
    step = _add_step(client, "FLASH_CHROMATOGRAPHY", category="PURIFICATION")

    response = _save(client, step, materials=[material])

    assert response.status_code == 400
    assert response.json["message"]


def test_database_refuses_a_material_with_two_links(
    app: Flask, client: FlaskClient, workup_setup: Dict[str, int]
) -> None:
    """The check constraint allows at most one link, even if the service is bypassed."""
    login(client)
    step = _add_step(client, "DRY_WITH_MATERIAL")

    with app.app_context():
        db.session.add(
            models.WorkupStepMaterial(
                step_id=step["id"],
                slot="drying_agent",
                sort_order=1,
                role=models.WorkupMaterialRole.DRYING_AGENT,
                compound_id=workup_setup["testanol_id"],
                solvent_name="Workupane",
                name="Two links",
                unit="g",
            )
        )
        with pytest.raises(IntegrityError, match="ck_workup_step_material_one_source"):
            db.session.commit()
        db.session.rollback()


def test_novel_compound_from_a_workbook_the_user_is_not_in_is_rejected(
    client: FlaskClient, workup_setup: Dict[str, int]
) -> None:
    """A novel compound id from someone else's workbook can't be linked."""
    login(client)
    step = _add_step(client, "DRY_WITH_MATERIAL")

    response = _save(
        client,
        step,
        materials=[
            {"slot": "drying_agent", "novelCompoundId": workup_setup["secretol_id"]}
        ],
    )

    assert response.status_code == 400


def test_single_slots_take_one_material_but_multiple_slots_take_mixtures(
    client: FlaskClient,
) -> None:
    """An eluent can be a mixture; a stationary phase is one material."""
    login(client)
    step = _add_step(client, "FLASH_CHROMATOGRAPHY", category="PURIFICATION")

    mixture = _save(
        client,
        step,
        materials=[
            {"slot": "eluent", "name": "Hexane", "quantity": 800},
            {"slot": "eluent", "name": "EtOAc", "quantity": 200},
        ],
    )
    two_silicas = _save(
        client,
        step,
        materials=[
            {"slot": "stationary_phase", "name": "Silica"},
            {"slot": "stationary_phase", "name": "Alumina"},
        ],
    )

    assert mixture.status_code == 200, mixture.json
    assert [m["name"] for m in mixture.json["materials"]] == ["Hexane", "EtOAc"]
    assert two_silicas.status_code == 400


def test_eluent_percentages_are_saved_with_each_eluent(client: FlaskClient) -> None:
    """Each eluent keeps its own percentage of the mixture."""
    login(client)
    step = _add_step(client, "FLASH_CHROMATOGRAPHY", category="PURIFICATION")

    response = _save(
        client,
        step,
        materials=[
            {"slot": "eluent", "name": "Hexane", "proportion": "90"},
            {"slot": "eluent", "name": "EtOAc", "proportion": 10},
        ],
    )

    assert response.status_code == 200, response.json
    assert [(m["name"], m["proportion"]) for m in response.json["materials"]] == [
        ("Hexane", 90),
        ("EtOAc", 10),
    ]


def test_changing_only_a_percentage_is_saved(client: FlaskClient) -> None:
    """A percentage-only change counts as a change, so it isn't skipped."""
    login(client)
    step = _add_step(client, "FLASH_CHROMATOGRAPHY", category="PURIFICATION")
    eluents = [
        {"slot": "eluent", "name": "Hexane", "proportion": 90},
        {"slot": "eluent", "name": "EtOAc", "proportion": 10},
    ]
    _save(client, step, materials=eluents)

    eluents[0]["proportion"] = 80
    eluents[1]["proportion"] = 20
    _save(client, step, materials=eluents)

    saved = _steps(client, category="PURIFICATION")["steps"][0]["materials"]
    assert [m["proportion"] for m in saved] == [80, 20]


@pytest.mark.parametrize(
    "proportions",
    [
        ["90", "20"],  # adds up to 110
        ["90", ""],  # only some given
        ["-10", "110"],
        ["ninety", "10"],
    ],
)
def test_invalid_eluent_percentages_are_rejected(
    client: FlaskClient, proportions: List[str]
) -> None:
    """Percentages must each be 0-100, all be given (or none) and add up to 100."""
    login(client)
    step = _add_step(client, "FLASH_CHROMATOGRAPHY", category="PURIFICATION")

    response = _save(
        client,
        step,
        materials=[
            {"slot": "eluent", "name": name, "proportion": proportion}
            for name, proportion in zip(("Hexane", "EtOAc"), proportions)
        ],
    )

    assert response.status_code == 400
    assert response.json["message"]


def test_percentages_are_ignored_for_slots_without_them(client: FlaskClient) -> None:
    """Only mixture slots (e.g. eluents) record a percentage."""
    login(client)
    step = _add_step(client, "FLASH_CHROMATOGRAPHY", category="PURIFICATION")

    response = _save(
        client,
        step,
        materials=[{"slot": "stationary_phase", "name": "Silica", "proportion": 50}],
    )

    assert response.status_code == 200, response.json
    assert response.json["materials"][0]["proportion"] is None


def test_fate_defaults_to_the_slot_default(client: FlaskClient) -> None:
    """The organic phase is kept by default, everything else is waste."""
    login(client)
    step = _add_step(client, "WASH")

    defaults = _save(
        client,
        step,
        materials=[
            {"slot": "organic_phase", "name": "EtOAc"},
            {"slot": "aqueous_phase", "name": "Brine"},
        ],
    )
    recovered = _save(
        client,
        step,
        materials=[{"slot": "organic_phase", "name": "EtOAc", "fate": "RECOVERED"}],
    )

    fates = {m["name"]: m["fate"] for m in defaults.json["materials"]}
    assert fates == {"EtOAc": "KEPT", "Brine": "WASTE"}
    assert recovered.json["materials"][0]["fate"] == "RECOVERED"


def test_total_quantity_counts_repeats_for_repeated_slots(
    app: Flask, client: FlaskClient
) -> None:
    """Three 20 mL washes are 60 mL of aqueous phase; the organic phase is used once."""
    login(client)
    step = _add_step(client, "WASH")
    _save(
        client,
        step,
        data={"repeats": 3},
        materials=[
            {"slot": "organic_phase", "name": "EtOAc", "quantity": 50},
            {"slot": "aqueous_phase", "name": "Brine", "quantity": 20},
        ],
    )

    with app.app_context():
        saved = db.session.get(models.WorkupStep, step["id"])
        totals = {
            m.slot: services.workup_steps.total_quantity(m) for m in saved.materials
        }
        saved.data = {}
        brine = next(m for m in saved.materials if m.slot == "aqueous_phase")
        without_repeats = services.workup_steps.total_quantity(brine)
        db.session.rollback()

    assert totals == {
        "organic_phase": (Decimal("50"), "mL"),
        "aqueous_phase": (Decimal("60"), "mL"),
    }
    assert without_repeats == (Decimal("20"), "mL")


def test_temperature_is_saved_and_validated(
    client: FlaskClient, audit_messages: List[Any]
) -> None:
    """Temperatures can be below zero, must be physical, and only where defined."""
    login(client)
    cooling = _add_step(client, "TEMPERATURE")
    wash = _add_step(client, "WASH")

    saved = _save(client, cooling, temperature="-78")
    too_cold = _save(client, cooling, temperature=-300)
    nonsense = _save(client, cooling, temperature="cold")
    not_defined = _save(client, wash, temperature=20)

    assert saved.status_code == 200, saved.json
    assert saved.json["temperature"] == -78
    assert too_cold.status_code == 400
    assert nonsense.status_code == 400
    assert not_defined.status_code == 400
    edit = next(m for m in audit_messages if m.field_name == "Workup Step Edited")
    assert edit.change_details["changes"]["temperature"] == {
        "old_value": None,
        "new_value": -78.0,
    }


def test_material_changes_are_audited(
    client: FlaskClient, audit_messages: List[Any]
) -> None:
    """The audit diff lists each slot's materials, and only when they change."""
    login(client)
    step = _add_step(client, "DRY_WITH_MATERIAL")
    materials = [{"slot": "drying_agent", "name": "MgSO4", "quantity": "2.5"}]

    _save(client, step, materials=materials)
    _save(client, step, materials=materials)

    edits = [m for m in audit_messages if m.field_name == "Workup Step Edited"]
    assert len(edits) == 1
    assert edits[0].change_details["changes"] == {
        "materials.drying_agent": {
            "old_value": None,
            "new_value": ["MgSO4 2.5 g (waste)"],
        }
    }


def test_solvent_list_has_chem21_flags_links_and_colours(
    client: FlaskClient, workup_setup: Dict[str, int]
) -> None:
    """The solvent picker gets the workbook's solvents and the user's colours."""
    login(client)

    response = client.get(f"/api/workups/{WORKGROUP}/{WORKBOOK}/solvents")

    assert response.status_code == 200
    by_name = {s["name"]: s for s in response.json["solvents"]}
    assert by_name["Testanol"]["flag"] == 2
    assert by_name["Testanol"]["compoundId"] == workup_setup["testanol_id"]
    assert by_name["Workupane"]["solventName"] == "Workupane"
    assert by_name["Workupane"]["compoundId"] is None
    assert response.json["colours"]["Recommended"]


def test_compound_lookup_finds_compounds_and_suggests(client: FlaskClient) -> None:
    """The compound picker uses the inventory's lookup."""
    login(client)
    url = f"/api/workups/{WORKGROUP}/{WORKBOOK}/compounds"

    found = client.get(url, query_string={"q": "Testanol"})
    partial = client.get(url, query_string={"q": "Test"})
    empty = client.get(url, query_string={"q": " "})

    assert found.status_code == 200
    assert found.json["status"] == "found"
    assert found.json["compound"]["name"] == "Testanol"
    assert partial.json["status"] == "not_found"
    assert "Testanol" in partial.json["suggestions"]
    assert empty.status_code == 400


def test_material_lookups_are_for_members_only(client: FlaskClient) -> None:
    """Non-members get 404s from the solvent list and compound lookup."""
    login(client, "wu_outsider", PASSWORD)

    solvents = client.get(f"/api/workups/{WORKGROUP}/{WORKBOOK}/solvents")
    compounds = client.get(
        f"/api/workups/{WORKGROUP}/{WORKBOOK}/compounds", query_string={"q": "Test"}
    )

    assert solvents.status_code == 404
    assert compounds.status_code == 404


# --- Templates ----------------------------------------------------------------


def test_save_as_template_keeps_step_types_not_values(client: FlaskClient) -> None:
    """A template stores the step types in order, with no values."""
    login(client)
    wash = _add_step(client, "WASH")
    _add_step(client, "DRY_WITH_MATERIAL")
    _save(client, wash, data={"repeats": 2})

    response = client.post(
        f"{_url()}/templates",
        json={"category": "WORKUP", "name": " Aqueous workup ", "description": ""},
    )
    listed = client.get(TEMPLATES, query_string={"category": "WORKUP"})

    assert response.status_code == 201, response.json
    assert response.json["name"] == "Aqueous workup"
    assert response.json["steps"] == ["WASH", "DRY_WITH_MATERIAL"]
    assert [t["name"] for t in listed.json] == ["Aqueous workup"]
    assert listed.json[0]["canDelete"] is True


def test_save_as_template_needs_a_name_and_steps(client: FlaskClient) -> None:
    """Templates can't be nameless or empty."""
    login(client)
    _add_step(client, "WASH")

    nameless = client.post(
        f"{_url()}/templates", json={"category": "WORKUP", "name": "  "}
    )
    empty = client.post(
        f"{_url()}/templates", json={"category": "PURIFICATION", "name": "Columns"}
    )

    assert nameless.status_code == 400
    assert empty.status_code == 400


def test_applying_a_template_appends_empty_steps(
    client: FlaskClient, audit_messages: List[Any]
) -> None:
    """Starting from a template copies its steps; the copies have no values."""
    login(client)
    wash = _add_step(client, "WASH")
    _add_step(client, "CONCENTRATION")
    _save(client, wash, data={"repeats": 2})
    template = client.post(
        f"{_url()}/templates", json={"category": "WORKUP", "name": "Standard"}
    ).json
    _add_step(client, "FILTRATION", code="WU1-003")

    response = client.post(
        f"{_url('WU1-003')}/apply-template", json={"templateId": template["id"]}
    )

    assert response.status_code == 200, response.json
    steps = response.json["steps"]
    assert [s["stepType"] for s in steps] == ["FILTRATION", "WASH", "CONCENTRATION"]
    assert all(s["data"] == {} for s in steps)
    assert "Workup Template Applied" in _field_names(audit_messages)


def test_template_from_another_workbook_cannot_be_applied(
    app: Flask, client: FlaskClient, workup_setup: Dict[str, int]
) -> None:
    """Templates are shared within one workbook only."""
    with app.app_context():
        other = models.WorkupTemplate(
            name="Elsewhere",
            workbook_id=workup_setup["other_workbook_id"],
            category=models.WorkupCategory.WORKUP,
            steps=[
                models.WorkupTemplateStep(
                    sort_order=1, step_type=models.WorkupStepType.WASH
                )
            ],
        )
        db.session.add(other)
        db.session.commit()
        other_id = other.id
    login(client)

    foreign = client.post(f"{_url()}/apply-template", json={"templateId": other_id})
    missing = client.post(f"{_url()}/apply-template", json={"templateId": 999999})
    nonsense = client.post(f"{_url()}/apply-template", json={"templateId": "abc"})
    delete = client.delete(f"{TEMPLATES}/{other_id}")

    assert foreign.status_code == 404
    assert missing.status_code == 404
    assert nonsense.status_code == 404
    assert delete.status_code == 404


def test_only_the_creator_can_delete_a_template(client: FlaskClient) -> None:
    """Other members can use a template but not delete it."""
    login(client)
    _add_step(client, "WASH")
    template = client.post(
        f"{_url()}/templates", json={"category": "WORKUP", "name": "Mine"}
    ).json

    login(client, "wu_member", PASSWORD)
    listed = client.get(TEMPLATES)
    refused = client.delete(f"{TEMPLATES}/{template['id']}")
    login(client)
    deleted = client.delete(f"{TEMPLATES}/{template['id']}")

    assert listed.json[0]["canDelete"] is False
    assert refused.status_code == 403
    assert deleted.status_code == 200
    assert client.get(TEMPLATES).json == []


# --- Access -------------------------------------------------------------------


def test_non_member_gets_404(client: FlaskClient) -> None:
    """Someone outside the workbook can't tell the reaction or templates exist."""
    login(client, "wu_outsider", PASSWORD)

    assert client.get(_url(), query_string={"category": "WORKUP"}).status_code == 404
    assert (
        client.post(
            f"{_url()}/steps", json={"category": "WORKUP", "stepType": "WASH"}
        ).status_code
        == 404
    )
    assert client.get(TEMPLATES).status_code == 404


def test_member_can_read_but_not_edit(client: FlaskClient) -> None:
    """Workbook members who didn't create the reaction get a read-only view."""
    login(client)
    step = _add_step(client, "WASH")
    login(client, "wu_member", PASSWORD)

    view = _steps(client)
    add = client.post(
        f"{_url()}/steps", json={"category": "WORKUP", "stepType": "WASH"}
    )
    update = client.put(f"{_url()}/steps/{step['id']}", json={"data": {}})
    delete = client.delete(f"{_url()}/steps/{step['id']}")

    assert view["canEdit"] is False
    assert len(view["steps"]) == 1
    assert add.status_code == 403
    assert update.status_code == 403
    assert delete.status_code == 403


def test_locked_reaction_is_read_only(client: FlaskClient) -> None:
    """Once a reaction is locked, even its creator can't change the workup."""
    login(client)

    view = _steps(client, code="WU1-002")
    add = client.post(
        f"{_url('WU1-002')}/steps", json={"category": "WORKUP", "stepType": "WASH"}
    )

    assert view["canEdit"] is False
    assert add.status_code == 403


def test_writes_need_json(client: FlaskClient) -> None:
    """Form posts are refused, which protects the routes from cross-site forms."""
    login(client)

    response = client.post(
        f"{_url()}/steps", data={"category": "WORKUP", "stepType": "WASH"}
    )

    assert response.status_code == 415
