"""
Tests for the ChemInventory integration (/api/inventory/<workgroup>/cheminventory).

ChemInventory itself is never called: `requests.post` is replaced with a fake
that answers like the ChemInventory API, and audit messages are captured
rather than sent.

Sets up its own workgroups so it doesn't affect other test modules:
- ChemInventory-Workgroup: the main test user ("test_username") is principal
  investigator, "ci_senior" a senior researcher and "ci_member" a standard member.
- ChemInventory-Other: owned by "ci_outsider".
"""

from typing import Any, Dict, List, Optional

import pytest
import requests
from cryptography.fernet import Fernet
from flask import Flask
from flask.testing import FlaskClient
from sources import models, services
from sources.extensions import db
from sources.services import cheminventory
from tests.populate_database import add_test_user
from tests.utils import add_workgroup, login

WORKGROUP = "ChemInventory-Workgroup"
OTHER_WORKGROUP = "ChemInventory-Other"
BASE = f"/api/inventory/{WORKGROUP}"
CONNECTION = f"{BASE}/cheminventory/connection"
SEARCH = f"{BASE}/cheminventory/search"
ACETONE_SEARCH = {"type": "cas", "query": "67-64-1"}

PASSWORD = "ci_pw"
ENCRYPTION_KEY = Fernet.generate_key().decode()
GOOD_TOKEN = "good-token-abcd1234"
ADMIN_TOKEN = "admin-token-9999"
BAD_TOKEN = "bad-token-0000"

# ChemInventory account types: 2 = Group Administrator, 6 = Read-only
ACCOUNT_TYPES = {GOOD_TOKEN: "6", ADMIN_TOKEN: "2"}

INVENTORY_ID = 9
CONTAINERS = [
    {
        "id": 101,
        "name": "Acetone",
        "location": 2,
        "locationtext": "Lab 1 > Cabinet A",
        "size": "500",
        "unit": "mL",
        "substance": 4454,  # a substance ID, not a name
        "cas": "67-64-1",
        "barcode": "CI-101",
        "supplier": "Sigma",
        "dateacquired": "2026-01-15",
        "comments": "private note",
    },
    {
        "id": 102,
        "name": "Sodium chloride",
        "location": 2,
        "locationtext": "Lab 1 > Cabinet A",
        "size": 1,
        "unit": "kg",
        "substance": 5234,
        "cas": "7647-14-5",
        "barcode": "CI-102",
        "supplier": "Fisher",
        "dateacquired": "15/01/2026",  # not ISO, so dropped
    },
    {
        "id": 103,
        "name": "Acetone (2)",
        "location": 3,
        "locationtext": "Lab 1 > Fridge",
        "size": "2.5",
        "unit": "L",
        "substance": 4454,  # a substance ID, not a name
        "cas": "67-64-1",
        "barcode": "CI-103",
        "supplier": "Sigma",
        "dateacquired": "",
    },
]


def _matches(container: Dict[str, Any], search_type: str, contents: str) -> bool:
    """Whether a fake container matches a search, roughly as ChemInventory does."""
    if search_type == "cas":
        return container["cas"] == contents
    return contents.lower() in container["name"].lower()


class FakeResponse:
    """Stands in for a requests.Response with a JSON body."""

    def __init__(self, payload: Any, status_code: int = 200) -> None:
        """
        Args:
            payload: The JSON body.
            status_code: The HTTP status.
        """
        self._payload = payload
        self.status_code = status_code
        self.ok = status_code < 400

    def json(self) -> Any:
        """Returns the JSON body."""
        return self._payload


class FakeChemInventory:
    """Answers ChemInventory API calls, and records them."""

    def __init__(self) -> None:
        self.calls: List[str] = []
        self.bodies: List[Dict[str, Any]] = []
        self.fail_with: Optional[Exception] = None
        self.respond_with: Optional[FakeResponse] = None
        self.containers: List[Dict[str, Any]] = CONTAINERS

    def post(
        self, url: str, json: Dict[str, Any], timeout: float, allow_redirects: bool
    ) -> FakeResponse:
        """Replaces requests.post."""
        assert allow_redirects is False, "redirects could forward the token"
        path = url.split("/api", 1)[1]
        self.calls.append(path)
        self.bodies.append({k: v for k, v in json.items() if k != "authtoken"})
        if self.fail_with is not None:
            raise self.fail_with
        if self.respond_with is not None:
            return self.respond_with
        token = json.get("authtoken")
        if token not in ACCOUNT_TYPES:
            return FakeResponse({"status": "error", "data": "Invalid token"})

        if path == "/general/getdetails":
            data: Any = {
                "user": {
                    "firstname": "Lab",
                    "lastname": "Manager",
                    "accounttype": ACCOUNT_TYPES[token],
                    "siteadmin": "0",
                    "enterpriseadmin": 0,
                },
                "inventory": {"id": INVENTORY_ID, "name": "Test Lab Inventory"},
            }
        elif path == "/search/execute":
            assert json["inventory"] == INVENTORY_ID
            data = {
                "orders": [],
                "containers": [
                    c
                    for c in self.containers
                    if _matches(c, json["type"], json["contents"])
                ],
            }
        else:
            raise AssertionError(f"Unexpected ChemInventory call: {path}")
        return FakeResponse({"status": "success", "data": data})


@pytest.fixture(scope="module")
def cheminventory_setup(app: Flask) -> Dict[str, int]:
    """
    Creates the users and workgroups for these tests.

    Returns:
        The main workgroup's id.
    """
    with app.app_context():
        pi = services.person.from_email("test_user@test.com")
        senior = add_test_user(
            "ci_senior", "ci_senior@test.com", "CI Senior", PASSWORD, True
        )
        member = add_test_user(
            "ci_member", "ci_member@test.com", "CI Member", PASSWORD, True
        )
        outsider = add_test_user(
            "ci_outsider", "ci_outsider@test.com", "CI Outsider", PASSWORD, True
        )
        workgroup = add_workgroup(
            WORKGROUP,
            principal_investigator=[pi],
            senior_researcher=[senior],
            standard_member=[member],
        )
        add_workgroup(OTHER_WORKGROUP, principal_investigator=[outsider])
        return {"workgroup_id": workgroup.id}


@pytest.fixture(autouse=True)
def fake_cheminventory(
    app: Flask, cheminventory_setup, monkeypatch: pytest.MonkeyPatch
) -> FakeChemInventory:
    """
    Switches the integration on, replaces ChemInventory with a fake, and
    starts each test with no connection.
    """
    monkeypatch.setitem(app.config, "CHEMINVENTORY_ENCRYPTION_KEY", ENCRYPTION_KEY)
    fake = FakeChemInventory()
    monkeypatch.setattr(cheminventory.requests, "post", fake.post)

    with app.app_context():
        db.session.query(models.ChemInventoryConnection).delete()
        db.session.commit()
    return fake


@pytest.fixture(autouse=True)
def audit_messages(monkeypatch: pytest.MonkeyPatch) -> List[Any]:
    """Captures audit messages instead of sending them."""
    sent: List[Any] = []
    monkeypatch.setattr(services.integration_history, "send_message", sent.append)
    return sent


def _connect(client: FlaskClient, token: Optional[str] = GOOD_TOKEN):
    """Saves a token for the workgroup and returns the response."""
    return client.put(CONNECTION, json={"token": token})


def _search(client: FlaskClient, body: Optional[Dict[str, Any]] = None):
    """Searches the workgroup's ChemInventory and returns the response."""
    return client.post(SEARCH, json=ACETONE_SEARCH if body is None else body)


def _stored_token(app: Flask) -> str:
    """The workgroup's token as stored in the database."""
    with app.app_context():
        return db.session.query(models.ChemInventoryConnection).one().encrypted_token


# --- Token encryption --------------------------------------------------------


def test_encrypt_round_trip(app: Flask):
    """Tokens are stored encrypted and can be decrypted with the same key."""
    with app.app_context():
        encrypted = cheminventory.encrypt_token(1, GOOD_TOKEN)
        assert GOOD_TOKEN not in encrypted
        assert cheminventory.decrypt_token(1, encrypted) == GOOD_TOKEN


def test_token_is_bound_to_its_workgroup(app: Flask):
    """A token encrypted for one workgroup can't be used by another."""
    with app.app_context():
        encrypted = cheminventory.encrypt_token(1, GOOD_TOKEN)
        with pytest.raises(cheminventory.ChemInventoryError):
            cheminventory.decrypt_token(2, encrypted)


def test_decrypt_with_removed_key_fails_cleanly(
    app: Flask, monkeypatch: pytest.MonkeyPatch
):
    """A token stored under a key that's been removed gives a readable error."""
    with app.app_context():
        encrypted = cheminventory.encrypt_token(1, GOOD_TOKEN)
        monkeypatch.setitem(
            app.config, "CHEMINVENTORY_ENCRYPTION_KEY", Fernet.generate_key().decode()
        )
        with pytest.raises(cheminventory.ChemInventoryError):
            cheminventory.decrypt_token(1, encrypted)


def test_key_rotation(app: Flask, client: FlaskClient, monkeypatch: pytest.MonkeyPatch):
    """
    With "new,old" keys, a token stored under the old key still works, and is
    re-encrypted with the new key when it's next used.
    """
    login(client)
    _connect(client)

    new_key = Fernet.generate_key().decode()
    monkeypatch.setitem(
        app.config,
        "CHEMINVENTORY_ENCRYPTION_KEY",
        f"{new_key},{ENCRYPTION_KEY}",  # noqa: E231
    )
    assert _search(client).status_code == 200

    monkeypatch.setitem(app.config, "CHEMINVENTORY_ENCRYPTION_KEY", new_key)
    assert _search(client).status_code == 200


def test_client_repr_hides_token():
    """The token never appears if the client is logged or printed."""
    client = cheminventory._Client(
        base_url="https://example.test", deadline=0, token=GOOD_TOKEN
    )
    assert GOOD_TOKEN not in repr(client)


# --- Configuration -----------------------------------------------------------


def test_switched_off_without_encryption_key(
    app: Flask, client: FlaskClient, monkeypatch: pytest.MonkeyPatch
):
    """Without an encryption key the integration says it's unavailable."""
    monkeypatch.setitem(app.config, "CHEMINVENTORY_ENCRYPTION_KEY", "")
    login(client)
    assert client.get(CONNECTION).json["available"] is False
    assert _search(client).status_code == 503
    assert _connect(client).status_code == 503


def test_switched_off_with_invalid_key(
    app: Flask, client: FlaskClient, monkeypatch: pytest.MonkeyPatch
):
    """An invalid encryption key switches the integration off."""
    monkeypatch.setitem(app.config, "CHEMINVENTORY_ENCRYPTION_KEY", "not-a-key")
    login(client)
    assert client.get(CONNECTION).json["available"] is False


def test_refuses_non_https_url(
    app: Flask,
    client: FlaskClient,
    fake_cheminventory: FakeChemInventory,
    monkeypatch: pytest.MonkeyPatch,
):
    """Tokens are never sent to a plain http:// URL."""
    monkeypatch.setitem(
        app.config, "CHEMINVENTORY_API_URL", "http://app.cheminventory.net/api"
    )
    login(client)
    assert client.get(CONNECTION).json["available"] is False
    assert _connect(client).status_code == 503
    assert fake_cheminventory.calls == []


# --- Access ------------------------------------------------------------------


def test_non_member_gets_404(client: FlaskClient):
    """Non-members can't see the connection or search."""
    login(client, "ci_outsider", PASSWORD)
    assert client.get(CONNECTION).status_code == 404
    assert _search(client).status_code == 404
    assert _connect(client).status_code == 404


def test_standard_member_can_read_but_not_manage(client: FlaskClient):
    """Standard members can view the status but not set or remove the token."""
    login(client, "ci_member", PASSWORD)
    response = client.get(CONNECTION)
    assert response.status_code == 200
    assert response.json["canManage"] is False

    assert _connect(client).status_code == 403
    assert client.delete(CONNECTION).status_code == 403


def test_senior_researcher_can_connect(client: FlaskClient):
    """Senior researchers can set the token."""
    login(client, "ci_senior", PASSWORD)
    assert _connect(client).status_code == 200


def test_writes_require_json(client: FlaskClient):
    """Form posts (as a cross-site form would send) are refused."""
    login(client)
    assert client.put(CONNECTION, data={"token": GOOD_TOKEN}).status_code == 415
    assert client.post(SEARCH, data=ACETONE_SEARCH).status_code == 415


def test_search_is_not_a_get(client: FlaskClient):
    """Searches can't be run with a GET, e.g. from another site's image tag."""
    login(client)
    assert client.get(SEARCH).status_code == 405


# --- Connection --------------------------------------------------------------


def test_status_when_not_connected(client: FlaskClient):
    """The status says the integration is available but not configured."""
    login(client)
    response = client.get(CONNECTION)
    assert response.json == {"available": True, "configured": False, "canManage": True}


def test_connect_stores_encrypted_token(
    app: Flask, client: FlaskClient, cheminventory_setup, audit_messages
):
    """
    Connecting checks the token, stores it encrypted, records the change, and
    never returns the token.
    """
    login(client)
    response = _connect(client)
    assert response.status_code == 200, response.json
    assert GOOD_TOKEN not in response.get_data(as_text=True)
    assert response.json["configured"] is True
    assert response.json["tokenHint"] == GOOD_TOKEN[-4:]
    assert response.json["cheminventoryUser"] == "Lab Manager"
    assert response.json["inventoryName"] == "Test Lab Inventory"
    assert response.json["connectedBy"]

    status = client.get(CONNECTION)
    assert GOOD_TOKEN not in status.get_data(as_text=True)

    stored = _stored_token(app)
    assert GOOD_TOKEN not in stored
    with app.app_context():
        workgroup_id = cheminventory_setup["workgroup_id"]
        assert cheminventory.decrypt_token(workgroup_id, stored) == GOOD_TOKEN

    assert [m.action for m in audit_messages] == ["connected"]
    assert audit_messages[0].token_hint == GOOD_TOKEN[-4:]
    assert GOOD_TOKEN not in audit_messages[0].serialise()


def test_replacing_token_is_audited(client: FlaskClient, audit_messages):
    """Replacing and removing the token are both recorded."""
    login(client)
    _connect(client)
    _connect(client)
    client.delete(CONNECTION)
    assert [m.action for m in audit_messages] == [
        "connected",
        "token_replaced",
        "disconnected",
    ]


def test_rejected_token_is_not_stored(app: Flask, client: FlaskClient, audit_messages):
    """A token ChemInventory refuses gives a 400 and nothing is saved."""
    login(client)
    response = _connect(client, BAD_TOKEN)
    assert response.status_code == 400
    assert "didn't accept" in response.json["message"]
    assert BAD_TOKEN not in response.get_data(as_text=True)
    with app.app_context():
        assert db.session.query(models.ChemInventoryConnection).count() == 0
    assert audit_messages == []


def test_admin_token_is_accepted(
    app: Flask, client: FlaskClient, caplog: pytest.LogCaptureFixture
):
    """
    Tokens for ChemInventory administrators are accepted, with a warning
    logged, as a read-only user's token is recommended.
    """
    login(client)
    response = _connect(client, ADMIN_TOKEN)
    assert response.status_code == 200, response.json
    with app.app_context():
        assert db.session.query(models.ChemInventoryConnection).count() == 1
    assert "administrator" in caplog.text
    assert ADMIN_TOKEN not in caplog.text


@pytest.mark.parametrize("token", ["", "   ", None, "has space", "x" * 600])
def test_malformed_token_is_refused(
    client: FlaskClient, fake_cheminventory: FakeChemInventory, token
):
    """Missing or malformed tokens are refused without calling ChemInventory."""
    login(client)
    assert _connect(client, token).status_code == 400
    assert fake_cheminventory.calls == []


def test_connect_when_unreachable(
    client: FlaskClient, fake_cheminventory: FakeChemInventory
):
    """If ChemInventory can't be reached, saving fails with a 502."""
    login(client)
    fake_cheminventory.fail_with = requests.ConnectionError("down")
    response = _connect(client)
    assert response.status_code == 502
    assert "couldn't be reached" in response.json["message"]


def test_redirects_are_not_followed(
    client: FlaskClient, fake_cheminventory: FakeChemInventory
):
    """A redirect is treated as an error, so the token is never forwarded."""
    login(client)
    fake_cheminventory.respond_with = FakeResponse({"status": "success"}, 307)
    assert _connect(client).status_code == 502


def test_disconnect(client: FlaskClient):
    """Disconnecting deletes the stored token."""
    login(client)
    _connect(client)
    assert client.delete(CONNECTION).status_code == 200
    assert client.get(CONNECTION).json["configured"] is False
    assert client.delete(CONNECTION).status_code == 409


# --- Search ------------------------------------------------------------------


def test_search_when_not_connected(client: FlaskClient):
    """Searching before connecting gives a 409 with a message."""
    login(client)
    response = _search(client)
    assert response.status_code == 409
    assert "hasn't connected" in response.json["message"]


def test_search_by_cas_groups_by_location(
    client: FlaskClient, fake_cheminventory: FakeChemInventory
):
    """
    A CAS search returns the matching containers grouped into one entry per
    location, sorted by name. Comments aren't shown, so aren't sent.
    """
    login(client)
    _connect(client)
    login(client, "ci_member", PASSWORD)

    response = _search(client)
    assert response.status_code == 200, response.json
    data = response.json
    body = response.get_data(as_text=True)
    assert GOOD_TOKEN not in body
    assert "private note" not in body
    assert data["type"] == "cas"
    assert data["query"] == "67-64-1"
    assert data["containerCount"] == 2
    assert data["truncated"] is False

    assert [loc["name"] for loc in data["locations"]] == [
        "Lab 1 > Cabinet A",
        "Lab 1 > Fridge",
    ]
    cabinet, fridge = data["locations"]
    assert cabinet["id"] == 2
    assert cabinet["path"] == ["Lab 1", "Cabinet A"]
    assert [c["id"] for c in cabinet["containers"]] == [101]
    container = cabinet["containers"][0]
    assert container["name"] == "Acetone"
    assert container["cas"] == "67-64-1"
    assert container["quantity"] == 500.0
    assert container["unit"] == "mL"
    assert container["barcode"] == "CI-101"
    assert container["dateReceived"] == "2026-01-15"
    assert fridge["containers"][0]["quantity"] == 2.5
    assert fridge["containers"][0]["dateReceived"] is None

    search_body = fake_cheminventory.bodies[-1]
    assert search_body == {
        "inventory": INVENTORY_ID,
        "type": "cas",
        "contents": "67-64-1",
    }


def test_search_by_name(client: FlaskClient):
    """A name search is passed to ChemInventory, trimmed."""
    login(client)
    _connect(client)
    response = _search(client, {"type": "name", "query": "  sodium  "})
    assert response.status_code == 200, response.json
    assert response.json["query"] == "sodium"
    [location] = response.json["locations"]
    assert location["containers"][0]["name"] == "Sodium chloride"
    assert location["containers"][0]["dateReceived"] is None  # not ISO, so dropped


def test_search_with_no_matches(client: FlaskClient):
    """A search that matches nothing returns no locations."""
    login(client)
    _connect(client)
    response = _search(client, {"type": "cas", "query": "50-00-0"})
    assert response.status_code == 200
    assert response.json["locations"] == []
    assert response.json["containerCount"] == 0


def test_search_asks_for_inventory_id_then_searches(
    client: FlaskClient, fake_cheminventory: FakeChemInventory
):
    """Each search asks which inventory the token belongs to, then searches it."""
    login(client)
    _connect(client)
    calls_before = len(fake_cheminventory.calls)
    _search(client)
    assert fake_cheminventory.calls[calls_before:] == [
        "/general/getdetails",
        "/search/execute",
    ]


def test_search_results_are_capped(
    client: FlaskClient,
    fake_cheminventory: FakeChemInventory,
    monkeypatch: pytest.MonkeyPatch,
):
    """Only the first MAX_SEARCH_RESULTS containers are kept, and that's flagged."""
    monkeypatch.setattr(cheminventory, "MAX_SEARCH_RESULTS", 1)
    login(client)
    _connect(client)
    response = _search(client)
    assert response.json["containerCount"] == 1
    assert response.json["truncated"] is True


@pytest.mark.parametrize(
    "body",
    [
        {},
        {"type": "barcode", "query": "CI-101"},
        {"type": "cas", "query": "acetone"},
        {"type": "cas", "query": "67-64-1; drop"},
        {"type": "name", "query": "a"},
        {"type": "name", "query": "x" * 101},
        {"type": "name", "query": 42},
    ],
)
def test_malformed_search_is_refused(
    client: FlaskClient, fake_cheminventory: FakeChemInventory, body
):
    """Malformed searches are refused without calling ChemInventory."""
    login(client)
    _connect(client)
    calls_before = len(fake_cheminventory.calls)
    assert _search(client, body).status_code == 400
    assert len(fake_cheminventory.calls) == calls_before


def test_search_time_limit(
    client: FlaskClient,
    fake_cheminventory: FakeChemInventory,
    monkeypatch: pytest.MonkeyPatch,
):
    """A search that runs past its time limit stops, with a 502."""
    login(client)
    _connect(client)
    monkeypatch.setattr(cheminventory, "SEARCH_TIME_LIMIT_SECONDS", 0)
    calls_before = len(fake_cheminventory.calls)

    response = _search(client)
    assert response.status_code == 502
    assert "too long" in response.json["message"]
    assert len(fake_cheminventory.calls) == calls_before


def test_unreachable_does_not_break_ai4green_inventory(
    client: FlaskClient, fake_cheminventory: FakeChemInventory
):
    """A ChemInventory outage gives a 502 there, and the AI4Green inventory still loads."""
    login(client)
    _connect(client)
    fake_cheminventory.fail_with = requests.ConnectionError("down")

    response = _search(client)
    assert response.status_code == 502
    assert "couldn't be reached" in response.json["message"]
    assert client.get(BASE).status_code == 200
