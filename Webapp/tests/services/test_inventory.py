"""
Tests for the inventory API (/api/inventory).

Sets up its own workgroups so it doesn't affect other test modules:
- Inventory-Workgroup: the main test user ("test_username") is principal
  investigator, "inv_senior" a senior researcher and "inv_member" a standard
  member. Inventory-Workbook belongs to it, with the main test user as member.
- Inventory-Other: owned by "inv_outsider", used to check one workgroup's ids
  can't reach another's data.
"""

from typing import Dict

import pytest
from flask import Flask
from flask.testing import FlaskClient
from sources import models
from sources.extensions import db
from tests.populate_database import add_test_user
from tests.utils import add_workgroup, login, main_test_person

WORKGROUP = "Inventory-Workgroup"
OTHER_WORKGROUP = "Inventory-Other"
BASE = f"/api/inventory/{WORKGROUP}"
OTHER_BASE = f"/api/inventory/{OTHER_WORKGROUP}"

# password shared by this module's users
PASSWORD = "inv_pw"

# seeded in tests/populate_database.py
COMPOUND_CAS = "123-45-5"
COMPOUND_NAME = "Testoic Acid"


@pytest.fixture(scope="module")
def inventory_setup(app: Flask) -> Dict[str, int]:
    """
    Creates the users, workgroups and workbook for these tests.

    Returns:
        Ids the tests need: the main workbook's id.
    """
    with app.app_context():
        pi = main_test_person()
        senior = add_test_user(
            "inv_senior", "inv_senior@test.com", "Senior Tester", PASSWORD, True
        )
        member = add_test_user(
            "inv_member", "inv_member@test.com", "Member Tester", PASSWORD, True
        )
        outsider = add_test_user(
            "inv_outsider", "inv_outsider@test.com", "Outside Tester", PASSWORD, True
        )

        workgroup = add_workgroup(
            WORKGROUP,
            principal_investigator=[pi],
            senior_researcher=[senior],
            standard_member=[member],
        )
        add_workgroup(OTHER_WORKGROUP, principal_investigator=[outsider])
        # no shared helper for workbooks yet; mirrors tests/populate_database.py
        workbook = models.WorkBook.create(
            name="Inventory-Workbook",
            group=workgroup.id,
            users=[pi],
            abbreviation="IW1",
        )
        db.session.commit()
        return {"workbook_id": workbook.id}


def _lookup(client: FlaskClient, query: str, base: str = BASE):
    """Looks up a compound in the workgroup's inventory and returns the response."""
    return client.get(f"{base}/compounds/lookup", query_string={"query": query})


def _add_location(client: FlaskClient, base: str, name: str) -> int:
    """Adds a location and returns its id."""
    response = client.post(f"{base}/locations", json={"name": name})
    assert response.status_code == 201, response.json
    return response.json["id"]


def _add_item(
    client: FlaskClient,
    base: str,
    location_id: int,
    cas: str = COMPOUND_CAS,
    **overrides,
) -> int:
    """Adds a container of a seeded compound (by default Testoic Acid) and returns its id."""
    lookup = _lookup(client, cas, base)
    compound = lookup.json["compound"]
    body = {
        "workgroupCompoundId": compound["workgroupCompoundId"],
        "compound": None
        if compound["workgroupCompoundId"]
        else {"source": compound["source"], "id": compound["id"]},
        "locationId": location_id,
        "quantity": 500,
        "unit": "mL",
        "supplier": "Test Supplies",
        "dateReceived": "2026-01-01",
        **overrides,
    }
    response = client.post(f"{base}/items", json=body)
    assert response.status_code == 201, response.json
    return response.json["id"]


# --- Access ------------------------------------------------------------------


def test_requires_login(client: FlaskClient, inventory_setup):
    """Anonymous users are redirected to log in rather than given data."""
    response = client.get(BASE)
    assert response.status_code == 302


def test_non_member_gets_404(client: FlaskClient, inventory_setup):
    """Non-members can't see the inventory, or tell the workgroup exists."""
    login(client, "inv_outsider", PASSWORD)
    assert client.get(BASE).status_code == 404
    assert client.get("/api/inventory/No-Such-Workgroup").status_code == 404


def test_standard_member_can_read_but_not_write(client: FlaskClient, inventory_setup):
    """Standard members get the inventory with canEdit false, and 403 on writes."""
    login(client, "inv_member", PASSWORD)
    response = client.get(BASE)
    assert response.status_code == 200
    assert response.json["permissions"] == {"canEdit": False}

    response = client.post(f"{BASE}/locations", json={"name": "Member shelf"})
    assert response.status_code == 403
    assert "senior researchers" in response.json["message"]


def test_senior_researcher_can_write(client: FlaskClient, inventory_setup):
    """Senior researchers can change the inventory."""
    login(client, "inv_senior", PASSWORD)
    assert client.get(BASE).json["permissions"] == {"canEdit": True}
    _add_location(client, BASE, "Senior shelf")


def test_writes_reject_form_posts(client: FlaskClient, inventory_setup):
    """Form-encoded bodies, which a cross-site form could send, are refused."""
    login(client)
    response = client.post(f"{BASE}/locations", data={"name": "Forged shelf"})
    assert response.status_code == 415


# --- Containers --------------------------------------------------------------


def test_add_adjust_move_and_dispose_container(client: FlaskClient, inventory_setup):
    """A container's full lifecycle, with each quantity change in the ledger."""
    login(client)
    cabinet = _add_location(client, BASE, "Lifecycle cabinet")
    fridge = _add_location(client, BASE, "Lifecycle fridge")
    item_id = _add_item(client, BASE, cabinet)

    response = client.post(
        f"{BASE}/items/{item_id}/adjust", json={"quantity": 350.5, "notes": "Used"}
    )
    assert response.status_code == 200
    response = client.post(f"{BASE}/items/{item_id}/move", json={"locationId": fridge})
    assert response.status_code == 200
    response = client.post(
        f"{BASE}/items/{item_id}/dispose", json={"notes": "Contaminated"}
    )
    assert response.status_code == 200

    inventory = client.get(BASE).json
    items = {item["id"]: item for wc in inventory["compounds"] for item in wc["items"]}
    item = items[item_id]
    assert item["isArchived"] is True
    assert item["quantity"] == 0
    assert item["location"]["id"] == fridge
    assert item["dateOpened"] is not None

    history = client.get(f"{BASE}/items/{item_id}/history").json
    assert [(e["reason"], e["change"]) for e in history] == [
        ("received", 500),
        ("adjustment", -149.5),
        ("disposed", -350.5),
    ]
    assert history[0]["user"]["name"] == "Gloria Testeban"

    # disposed containers can't be changed again
    response = client.post(f"{BASE}/items/{item_id}/adjust", json={"quantity": 1})
    assert response.status_code == 400


@pytest.mark.parametrize(
    "overrides",
    [
        {"quantity": -1},
        {"quantity": 0},
        {"quantity": "500"},
        {"quantity": True},
        {"quantity": 1e12},
        {"unit": "bucket"},
        {"expiryDate": "not a date"},
        {"lotNumber": "x" * 101},
    ],
)
def test_add_container_validates_input(
    client: FlaskClient, inventory_setup, overrides, request: pytest.FixtureRequest
):
    """Invalid fields are rejected with a 400 and a message."""
    login(client)
    # one location per case, named after the case (e.g. "overrides3")
    location_id = _add_location(client, BASE, f"Validation {request.node.callspec.id}")
    lookup = _lookup(client, COMPOUND_CAS)
    compound = lookup.json["compound"]
    body = {
        "compound": {"source": compound["source"], "id": compound["id"]},
        "locationId": location_id,
        "quantity": 10,
        "unit": "g",
        **overrides,
    }
    response = client.post(f"{BASE}/items", json=body)
    assert response.status_code == 400
    assert response.json["message"]


def test_archived_location_blocked_while_it_has_containers(
    client: FlaskClient, inventory_setup
):
    """A location with active containers can't be archived; an empty one can."""
    login(client)
    full = _add_location(client, BASE, "Full shelf")
    empty = _add_location(client, BASE, "Empty shelf")
    _add_item(client, BASE, full)

    response = client.patch(f"{BASE}/locations/{full}", json={"isArchived": True})
    assert response.status_code == 400
    response = client.patch(f"{BASE}/locations/{empty}", json={"isArchived": True})
    assert response.status_code == 200 and response.json["isArchived"] is True

    # archived locations can't receive containers
    lookup = _lookup(client, COMPOUND_CAS)
    response = client.post(
        f"{BASE}/items",
        json={
            "workgroupCompoundId": lookup.json["compound"]["workgroupCompoundId"],
            "locationId": empty,
            "quantity": 1,
            "unit": "g",
        },
    )
    assert response.status_code == 400


def test_barcodes_are_unique_among_active_containers(
    client: FlaskClient, inventory_setup
):
    """A barcode already on an active container is refused, until that
    container is disposed of."""
    login(client)
    shelf = _add_location(client, BASE, "Barcode shelf")
    item_id = _add_item(client, BASE, shelf, barcode="INV-BC-1")

    lookup = _lookup(client, COMPOUND_CAS)
    body = {
        "workgroupCompoundId": lookup.json["compound"]["workgroupCompoundId"],
        "locationId": shelf,
        "quantity": 1,
        "unit": "g",
        "barcode": " INV-BC-1 ",  # trimmed before it's compared
    }
    response = client.post(f"{BASE}/items", json=body)
    assert response.status_code == 400
    assert "INV-BC-1" in response.json["message"]

    client.post(f"{BASE}/items/{item_id}/dispose", json={})
    response = client.post(f"{BASE}/items", json=body)
    assert response.status_code == 201


def test_location_names_are_unique_ignoring_case(client: FlaskClient, inventory_setup):
    """A second location with the same name in a different case is refused."""
    login(client)
    _add_location(client, BASE, "Unique Cabinet")
    response = client.post(f"{BASE}/locations", json={"name": "unique cabinet"})
    assert response.status_code == 400


def test_ids_from_another_workgroup_are_not_found(client: FlaskClient, inventory_setup):
    """A member of one workgroup can't touch another workgroup's rows by id."""
    login(client, "inv_outsider", PASSWORD)
    other_location = _add_location(client, OTHER_BASE, "Outsider shelf")
    other_item = _add_item(client, OTHER_BASE, other_location)

    client.get("/auth/logout")
    login(client)
    assert client.get(f"{BASE}/items/{other_item}/history").status_code == 404
    response = client.post(f"{BASE}/items/{other_item}/adjust", json={"quantity": 1})
    assert response.status_code == 404
    response = client.put(
        f"{BASE}/locations/{other_location}", json={"name": "Taken over"}
    )
    assert response.status_code == 404

    # nor put its own container into the other workgroup's location
    own_location = _add_location(client, BASE, "Own shelf")
    own_item = _add_item(client, BASE, own_location)
    response = client.post(
        f"{BASE}/items/{own_item}/move", json={"locationId": other_location}
    )
    assert response.status_code == 404


# --- Lookup and novel compounds ----------------------------------------------


def test_lookup_by_cas_and_name(client: FlaskClient, inventory_setup):
    """CAS and exact-name lookups find the compound; unknown ones say so."""
    login(client)

    by_cas = _lookup(client, COMPOUND_CAS)
    assert by_cas.json["status"] == "found"
    assert by_cas.json["compound"]["name"] == COMPOUND_NAME
    assert by_cas.json["compound"]["hazards"] == ["H900", "H901"]
    assert by_cas.json["compound"]["hazardStatements"] == {
        "H900": "Caution testing",
        "H901": "May cause tests",
    }

    by_name = _lookup(client, COMPOUND_NAME.lower())
    assert by_name.json["compound"]["cas"] == COMPOUND_CAS

    unknown_cas = _lookup(client, "9999-99-9")
    assert unknown_cas.json == {"status": "cas_not_found", "cas": "9999-99-9"}

    partial = _lookup(client, "Testa")
    assert partial.json["status"] == "not_found"
    assert "Testamine" in partial.json["suggestions"]

    # % is matched literally, not as a SQL wildcard
    wildcard = _lookup(client, "%")
    assert wildcard.json == {"status": "not_found", "suggestions": []}

    assert client.get(f"{BASE}/compounds/lookup").status_code == 400


def test_add_novel_compound_and_stock_it(
    client: FlaskClient, inventory_setup: Dict[str, int]
):
    """A novel compound can be added to a workbook, found by CAS and stocked."""
    login(client)
    response = client.post(
        f"{BASE}/novel-compounds",
        json={
            "workbookId": inventory_setup["workbook_id"],
            "name": "Inventorium",
            "cas": "7777-77-7",
            "smiles": "CCCCCCO",
            "hazardCodes": "",
        },
    )
    assert response.status_code == 201, response.json
    compound = response.json
    assert compound["source"] == "novel"
    assert compound["workbook"] == "Inventory-Workbook"

    lookup = _lookup(client, "7777-77-7")
    assert lookup.json["compound"]["id"] == compound["id"]

    location_id = _add_location(client, BASE, "Novel shelf")
    response = client.post(
        f"{BASE}/items",
        json={
            "compound": {"source": "novel", "id": compound["id"]},
            "locationId": location_id,
            "quantity": 2,
            "unit": "g",
        },
    )
    assert response.status_code == 201

    names = [wc["compound"]["name"] for wc in client.get(BASE).json["compounds"]]
    assert "Inventorium" in names


def test_novel_compound_needs_a_workbook_the_user_belongs_to(
    client: FlaskClient, inventory_setup: Dict[str, int]
):
    """Senior researchers who aren't in the workbook can't add compounds to it."""
    login(client, "inv_senior", PASSWORD)
    response = client.post(
        f"{BASE}/novel-compounds",
        json={
            "workbookId": inventory_setup["workbook_id"],
            "name": "Sneakium",
            "cas": "8888-88-8",
        },
    )
    assert response.status_code == 404


def test_novel_compound_validation_errors_are_reported(
    client: FlaskClient, inventory_setup: Dict[str, int]
):
    """Validation from the reaction table's novel compound code is returned as a message."""
    login(client)
    response = client.post(
        f"{BASE}/novel-compounds",
        json={
            "workbookId": inventory_setup["workbook_id"],
            "name": "Bad smiles",
            "smiles": "not-a-smiles((",
        },
    )
    assert response.status_code == 400
    assert "SMILES" in response.json["message"]


# --- Structure search --------------------------------------------------------

# seeded in tests/populate_database.py: Testoic Acid (COMPOUND_CAS) is methane
# (C) and Testamine is ethane (CC)
ETHANE_CAS = "123-45-7"


def _structure_search(client: FlaskClient, smiles: str, mode: str, base: str = BASE):
    """Runs a structure search on the workgroup's inventory and returns the response."""
    return client.post(
        f"{base}/structure-search", json={"smiles": smiles, "mode": mode}
    )


def _workgroup_compound_id(client: FlaskClient, cas: str) -> int:
    """The workgroup's WorkgroupCompound id for a stocked compound."""
    return _lookup(client, cas).json["compound"]["workgroupCompoundId"]


def test_structure_search(client: FlaskClient, inventory_setup):
    """Exact searches find the same compound; substructure searches find compounds containing it."""
    login(client)
    location_id = _add_location(client, BASE, "Structure search shelf")
    _add_item(client, BASE, location_id)
    _add_item(client, BASE, location_id, cas=ETHANE_CAS)
    methane = _workgroup_compound_id(client, COMPOUND_CAS)
    ethane = _workgroup_compound_id(client, ETHANE_CAS)

    response = _structure_search(client, "C", "exact")
    assert response.status_code == 200
    ids = response.json["workgroupCompoundIds"]
    assert methane in ids and ethane not in ids

    ids = _structure_search(client, "CC", "substructure").json["workgroupCompoundIds"]
    assert ethane in ids and methane not in ids

    ids = _structure_search(client, "C", "substructure").json["workgroupCompoundIds"]
    assert methane in ids and ethane in ids


def test_standard_member_can_structure_search(client: FlaskClient, inventory_setup):
    """Structure search only reads, so any member can use it."""
    login(client, "inv_member", PASSWORD)
    assert _structure_search(client, "C", "exact").status_code == 200


def test_structure_search_non_member_gets_404(client: FlaskClient, inventory_setup):
    """Non-members can't search the inventory, or tell the workgroup exists."""
    login(client, "inv_outsider", PASSWORD)
    assert _structure_search(client, "C", "exact").status_code == 404


def test_structure_search_rejects_form_posts(client: FlaskClient, inventory_setup):
    """Only JSON bodies are accepted."""
    login(client)
    response = client.post(
        f"{BASE}/structure-search", data={"smiles": "C", "mode": "exact"}
    )
    assert response.status_code == 415


@pytest.mark.parametrize(
    "smiles, mode, message",
    [
        ("thisIsDefinitelyNotASmilesString", "exact", "Invalid structure"),
        ("", "substructure", "Draw a structure"),
        ("C", "similar", "exact or a substructure"),
        ("*C*", "exact", "Polymer structures"),
    ],
)
def test_structure_search_validation(
    client: FlaskClient, inventory_setup, smiles: str, mode: str, message: str
):
    """Structures that can't be searched are a 400 with a message to show."""
    login(client)
    response = _structure_search(client, smiles, mode)
    assert response.status_code == 400
    assert message in response.json["message"]
