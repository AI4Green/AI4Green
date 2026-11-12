"""
Tests for the reaction search API (/api/search).

Uses the reactions seeded in tests/populate_database.py, in Test-Workgroup /
Test-Workbook:
- TW1-001 "a test reaction": reactants C, CC; reagents CCC, CCCCC
  (Testide, Testol); solvent Testanol; products CP, CPP
- TW1-002 "a test polymer reaction": the same, with the polymer *C* in
  place of C

Other test modules can add reactions to the same workbook, so these tests
check which reactions are or aren't found rather than exact counts.
"""

from datetime import datetime, timedelta
from typing import Dict, Iterator, List

import pytest
import pytz
from flask import Flask
from flask.testing import FlaskClient
from sources import models
from sources.extensions import db
from tests.populate_database import add_test_user
from tests.utils import login

SEARCH_URL = "/api/search/reactions"
WORKGROUP = "Test-Workgroup"
WORKBOOK = "Test-Workbook"
REACTION = "TW1-001"
POLYMER_REACTION = "TW1-002"
OUTSIDER = "search_outsider"
OUTSIDER_PASSWORD = "search_pw"
# hexamine: in none of the seeded reactions
UNUSED_SMILES = "C1N2CN3CN1CN(C2)C3"


def _search(client: FlaskClient, **criteria) -> Dict:
    """Searches the Test-Workbook and returns the JSON response, checking it succeeded."""
    body = {"workgroup": WORKGROUP, "workbook": WORKBOOK, **criteria}
    response = client.post(SEARCH_URL, json=body)
    assert response.status_code == 200, response.json
    return response.json


def _structure(smiles: str, mode: str = "exact", roles: List[str] = None) -> Dict:
    """Builds the structure part of a search."""
    structure = {"smiles": smiles, "mode": mode}
    if roles is not None:
        structure["roles"] = roles
    return structure


def _found(response: Dict) -> Dict[str, Dict]:
    """The results keyed by reaction ID."""
    return {result["reactionId"]: result for result in response["results"]}


@pytest.fixture(scope="module")
def outsider(app: Flask) -> None:
    """Adds a verified user who isn't in Test-Workgroup."""
    with app.app_context():
        add_test_user(
            OUTSIDER, f"{OUTSIDER}@test.com", "Search Outsider", OUTSIDER_PASSWORD, True
        )
        db.session.commit()


@pytest.fixture
def reaction_with_testanol_solvent(app: Flask) -> Iterator[None]:
    """
    Points TW1-001's solvent at the Testanol compound's primary key for one
    test, then restores it. The seeded value (-5) isn't a valid key.
    """
    with app.app_context():
        reaction = _seeded_reaction(REACTION)
        original_solvent = list(reaction.solvent)
        testanol = db.session.query(models.Compound).filter_by(cas="123-46-0").one()
        reaction.solvent = [str(testanol.id)]
        db.session.commit()
    yield
    with app.app_context():
        reaction = _seeded_reaction(REACTION)
        reaction.solvent = original_solvent
        db.session.commit()


def _seeded_reaction(reaction_id: str) -> models.Reaction:
    """Gets a seeded reaction from Test-Workbook by its reaction ID."""
    return (
        db.session.query(models.Reaction)
        .join(models.WorkBook, models.Reaction.workbooks == models.WorkBook.id)
        .filter(
            models.WorkBook.name == WORKBOOK, models.Reaction.reaction_id == reaction_id
        )
        .one()
    )


# --- Access and validation ---------------------------------------------------


def test_requires_login(client: FlaskClient):
    """Anonymous users are redirected to log in rather than given results."""
    response = client.post(SEARCH_URL, json={"keyword": "test"})
    assert response.status_code == 302


def test_rejects_form_posts(client: FlaskClient):
    """Only JSON bodies are accepted, so cross-site forms can't post to it."""
    login(client)
    response = client.post(SEARCH_URL, data={"keyword": "test"})
    assert response.status_code == 415


def test_needs_a_search_term(client: FlaskClient):
    """A workgroup or workbook alone isn't a search."""
    login(client)
    response = client.post(SEARCH_URL, json={"workgroup": WORKGROUP})
    assert response.status_code == 400
    assert "keyword, reagent, date or structure" in response.json["message"]


@pytest.mark.parametrize(
    "body, message",
    [
        ({"workbook": WORKBOOK, "keyword": "test"}, "Choose a workgroup"),
        (
            {"keyword": "test", "dateFrom": "2026-02-01", "dateTo": "2026-01-01"},
            "start date must be on or before",
        ),
        ({"dateFrom": "yesterday"}, "must be a date"),
        ({"keyword": "x" * 201}, "at most 200 characters"),
        ({"structure": _structure("C", mode="similar")}, "exact or a substructure"),
        ({"structure": _structure("C", roles=[])}, "at least one place"),
        ({"structure": _structure("C", roles=["catalyst"])}, "reactants, reagents"),
        (
            {"structure": _structure("thisIsDefinitelyNotASmilesString")},
            "Invalid structure, please try again!",
        ),
        (
            {"structure": _structure("*C*", mode="substructure")},
            "not supported for polymers",
        ),
    ],
)
def test_invalid_requests_are_reported(client: FlaskClient, body: Dict, message: str):
    """Invalid input is a 400 with a message the page can show."""
    login(client)
    response = client.post(SEARCH_URL, json=body)
    assert response.status_code == 400
    assert message in response.json["message"]


def test_unknown_workgroup_or_workbook_is_not_found(client: FlaskClient):
    """Workgroups and workbooks the user can't see get a 404."""
    login(client)
    for body in (
        {"workgroup": "No-Such-Workgroup", "keyword": "test"},
        {"workgroup": WORKGROUP, "workbook": "No-Such-Workbook", "keyword": "test"},
    ):
        assert client.post(SEARCH_URL, json=body).status_code == 404


def test_non_member_cannot_search_the_workgroup(client: FlaskClient, outsider):
    """Non-members get a 404 for the workgroup, and never see its reactions."""
    login(client, OUTSIDER, OUTSIDER_PASSWORD)
    response = client.post(SEARCH_URL, json={"workgroup": WORKGROUP, "keyword": "test"})
    assert response.status_code == 404

    response = client.post(SEARCH_URL, json={"keyword": "test"})
    assert response.status_code == 200
    assert all(result["workgroup"] != WORKGROUP for result in response.json["results"])


def test_lists_the_users_workgroups_and_workbooks(client: FlaskClient):
    """The filter options are the user's workgroups and their workbooks."""
    login(client)
    response = client.get("/api/search/workgroups")
    assert response.status_code == 200
    workgroups = {wg["name"]: wg["workbooks"] for wg in response.json["workgroups"]}
    assert WORKBOOK in workgroups[WORKGROUP]


# --- Exact structure search (the cases from the old Jinja search) -------------


def test_exact_structure_search_finds_the_reaction(client: FlaskClient):
    """Methane is a reactant in TW1-001; the polymer reaction has *C* instead."""
    login(client)
    found = _found(_search(client, structure=_structure("C")))
    assert REACTION in found
    assert POLYMER_REACTION not in found

    result = found[REACTION]
    assert result["name"] == "a test reaction"
    assert result["workgroup"] == WORKGROUP
    assert result["workbook"] == WORKBOOK
    assert result["matchedRoles"] == ["reactant"]
    assert result["creatorUsername"] == "test_username"
    # the seeded reactions have no scheme image
    assert result["image"] is None


@pytest.mark.parametrize(
    "scope",
    [
        {"workgroup": WORKGROUP},  # every workbook in the workgroup
        {},  # every workgroup
    ],
)
def test_exact_structure_search_across_workbooks_and_workgroups(
    client: FlaskClient, scope: Dict
):
    """Leaving out the workbook, or the workgroup too, widens the search."""
    login(client)
    response = client.post(SEARCH_URL, json={**scope, "structure": _structure("C")})
    assert response.status_code == 200
    assert REACTION in _found(response.json)


def test_exact_structure_search_with_no_results(client: FlaskClient):
    """A structure in no reaction finds nothing."""
    login(client)
    response = _search(client, structure=_structure(UNUSED_SMILES))
    assert response == {"results": [], "truncated": False}


def test_exact_polymer_structure_search(client: FlaskClient):
    """Polymer repeat units are found in polymer reactions only."""
    login(client)
    found = _found(_search(client, structure=_structure("*C*")))
    assert POLYMER_REACTION in found
    assert REACTION not in found


# --- Substructure search and roles --------------------------------------------


def test_substructure_search(client: FlaskClient):
    """
    Ethane (CC) is in the reactant ethane AND the reagents propane and pentane
    of both reactions, but not in the one-carbon products.
    """
    login(client)
    found = _found(_search(client, structure=_structure("CC", mode="substructure")))
    assert REACTION in found
    assert POLYMER_REACTION in found
    assert found[REACTION]["matchedRoles"] == ["reactant", "reagent"]


def test_substructure_search_in_reagents_only(client: FlaskClient):
    """Butane is only contained in the reagent pentane (CCCCC)."""
    login(client)
    found = _found(
        _search(
            client,
            structure=_structure("CCCC", mode="substructure", roles=["reagent"]),
        )
    )
    assert found[REACTION]["matchedRoles"] == ["reagent"]


def test_roles_limit_where_the_structure_is_matched(client: FlaskClient):
    """Propane is a reagent, so searching only products doesn't find it."""
    login(client)
    found = _found(_search(client, structure=_structure("CCC", roles=["reagent"])))
    assert found[REACTION]["matchedRoles"] == ["reagent"]

    found = _found(_search(client, structure=_structure("CCC", roles=["product"])))
    assert REACTION not in found


def test_structure_search_in_solvents(
    client: FlaskClient, reaction_with_testanol_solvent
):
    """Solvents are stored as compound keys and resolved to SMILES for matching."""
    login(client)
    found = _found(_search(client, structure=_structure("CCCCC", roles=["solvent"])))
    assert found[REACTION]["matchedRoles"] == ["solvent"]

    # pentane is also a reagent, so searching everywhere matches both
    found = _found(_search(client, structure=_structure("CCCCC")))
    assert found[REACTION]["matchedRoles"] == ["reagent", "solvent"]


def test_invalid_solvent_keys_are_skipped(client: FlaskClient):
    """The seeded solvent key (-5) doesn't resolve, and doesn't break the search."""
    login(client)
    found = _found(_search(client, structure=_structure("CCCCC", roles=["solvent"])))
    assert REACTION not in found


# --- Text and date filters ---------------------------------------------------


def test_keyword_search(client: FlaskClient):
    """Keywords match the name, description or reaction ID, ignoring case."""
    login(client)
    found = _found(_search(client, keyword="POLYMER"))
    assert POLYMER_REACTION in found
    assert REACTION not in found

    assert REACTION in _found(_search(client, keyword="tw1-001"))
    assert REACTION in _found(_search(client, keyword="test family"))


def test_keyword_wildcards_are_matched_literally(client: FlaskClient):
    """A LIKE wildcard in the keyword doesn't match everything."""
    login(client)
    assert _search(client, keyword="%")["results"] == []


def test_reagent_and_solvent_name_search(client: FlaskClient):
    """Reagent and solvent names match case-insensitively."""
    login(client)
    assert REACTION in _found(_search(client, reagent="testide"))  # reagent
    assert REACTION in _found(_search(client, reagent="TESTANOL"))  # solvent
    assert REACTION not in _found(_search(client, reagent="no such reagent"))


def test_date_range_search(client: FlaskClient):
    """The creation date range is inclusive."""
    login(client)
    today = datetime.now(pytz.timezone("Europe/London")).date()
    found = _found(
        _search(
            client,
            dateFrom=(today - timedelta(days=1)).isoformat(),
            dateTo=(today + timedelta(days=1)).isoformat(),
        )
    )
    assert REACTION in found

    found = _found(_search(client, dateFrom=(today + timedelta(days=2)).isoformat()))
    assert REACTION not in found


def test_filters_combine(client: FlaskClient):
    """Every given filter must match."""
    login(client)
    found = _found(_search(client, keyword="polymer", structure=_structure("CC")))
    assert POLYMER_REACTION in found
    assert REACTION not in found
