"""
Tests for the inventory CSV import (/api/inventory/<workgroup>/import).

Sets up its own workgroup so it doesn't affect other test modules:
Import-Workgroup, where the main test user is principal investigator and
"imp_member" a standard member, with Import-Workbook for novel compounds.
"imp_outsider" isn't a member.

libmagic is mocked to report each file's MIME type, so the tests don't depend
on the libmagic version installed, or on it being installed at all.
"""

import io
from typing import Dict, Optional

import pytest
import pytest_mock
from flask import Flask
from flask.testing import FlaskClient
from sources import models
from sources.extensions import db
from sources.services.inventory_import import MAX_ROWS
from tests.populate_database import add_test_user
from tests.utils import add_workgroup, login, main_test_person

WORKGROUP = "Import-Workgroup"
BASE = f"/api/inventory/{WORKGROUP}"
PASSWORD = "imp_pw"
HEADERS = {"X-AI4Green-Request": "inventory-import"}

# seeded in tests/populate_database.py
COMPOUND_CAS = "123-45-5"
COMPOUND_NAME = "Testoic Acid"

HEADER_ROW = (
    "CAS number,Name,SMILES,Hazard codes,Location,Location description,"
    "Quantity,Unit\n"
)
# structures no seeded compound has
HEXANEDIOL = "OCCCCCCO"
BUTANEDIOL = "OCCCCO"


def _row(
    cas: str,
    name: str,
    location: str,
    quantity: str,
    unit: str,
    smiles: str = "",
    hazard_codes: str = "",
    location_description: str = "",
) -> str:
    """One CSV line in HEADER_ROW's column order."""
    cells = [
        cas,
        name,
        smiles,
        hazard_codes,
        location,
        location_description,
        quantity,
        unit,
    ]
    return ",".join(cells) + "\n"


# a row of the seeded compound in "Shelf"
KNOWN_ROW = _row(COMPOUND_CAS, COMPOUND_NAME, "Shelf", "1", "g")


@pytest.fixture(scope="module")
def import_setup(app: Flask) -> Dict[str, int]:
    """
    Creates the users, workgroup and workbook for these tests.

    Returns:
        Ids the tests need: the workbook's id.
    """
    with app.app_context():
        pi = main_test_person()
        member = add_test_user(
            "imp_member", "imp_member@test.com", "Import Member", PASSWORD, True
        )
        add_test_user(
            "imp_outsider", "imp_outsider@test.com", "Import Outsider", PASSWORD, True
        )
        workgroup = add_workgroup(
            WORKGROUP, principal_investigator=[pi], standard_member=[member]
        )
        # no shared helper for workbooks yet; mirrors tests/populate_database.py
        workbook = models.WorkBook.create(
            name="Import-Workbook",
            group=workgroup.id,
            users=[pi],
            abbreviation="IMP",
        )
        db.session.commit()
        return {"workbook_id": workbook.id}


@pytest.fixture
def sniffed_type(mocker: pytest_mock.MockerFixture):
    """Makes libmagic report text/csv; tests can change `return_value`."""
    return mocker.patch(
        "sources.services.file_attachments.magic.Magic.from_buffer",
        return_value="text/csv",
    )


def _post(
    client: FlaskClient,
    content: bytes,
    preview: bool = False,
    filename: str = "stock.csv",
    workbook_id: Optional[int] = None,
    headers: Optional[Dict[str, str]] = None,
):
    """Uploads a file to the import or preview route and returns the response."""
    data = {"file": (io.BytesIO(content), filename)}
    if workbook_id is not None:
        data["workbookId"] = str(workbook_id)
    return client.post(
        f"{BASE}/import/preview" if preview else f"{BASE}/import",
        data=data,
        content_type="multipart/form-data",
        headers=HEADERS if headers is None else headers,
    )


def _inventory(client: FlaskClient) -> Dict:
    """The workgroup's inventory, as the page gets it."""
    return client.get(BASE).json


def _containers_of(client: FlaskClient, name: str) -> list:
    """The containers the workgroup holds of a compound."""
    for compound in _inventory(client)["compounds"]:
        if compound["compound"]["name"] == name:
            return compound["items"]
    return []


# --- Access and CSRF ---------------------------------------------------------


def test_non_member_gets_404(client: FlaskClient, import_setup, sniffed_type):
    """Non-members can't import, or tell the workgroup exists."""
    login(client, "imp_outsider", PASSWORD)
    assert _post(client, HEADER_ROW.encode()).status_code == 404


def test_standard_member_gets_403(client: FlaskClient, import_setup, sniffed_type):
    """Standard members can't import."""
    login(client, "imp_member", PASSWORD)
    assert _post(client, HEADER_ROW.encode()).status_code == 403


def test_missing_csrf_header_is_refused(
    client: FlaskClient, import_setup, sniffed_type
):
    """A multipart post without the custom header, as a forged form would send,
    is refused before the file is read."""
    login(client)
    content = (HEADER_ROW + KNOWN_ROW).encode()
    response = _post(client, content, headers={})
    assert response.status_code == 400
    assert _containers_of(client, COMPOUND_NAME) == []


# --- File validation ---------------------------------------------------------


@pytest.mark.parametrize(
    "content, filename, mime_type",
    [
        (b"\x89PNG\r\n\x1a\n", "stock.csv", "image/png"),  # a PNG renamed .csv
        (b"PK\x03\x04", "stock.csv", "application/zip"),
        (b"a,b\n", "stock.xlsx", "text/csv"),
        (b"a,b\n", "stock.txt", "text/plain"),
        (b"a" * (1024 * 1024 + 1), "stock.csv", "text/csv"),
    ],
)
def test_invalid_files_are_refused(
    client: FlaskClient,
    import_setup,
    sniffed_type,
    content: bytes,
    filename: str,
    mime_type: str,
):
    """Only CSV files of up to 1 MB are read, judged by the sniffed type."""
    login(client)
    sniffed_type.return_value = mime_type
    response = _post(client, content, preview=True, filename=filename)
    assert response.status_code == 400
    assert "CSV file" in response.json["message"]


def test_non_utf8_file_is_refused(client: FlaskClient, import_setup, sniffed_type):
    """Files that aren't UTF-8 are refused rather than read with garbled text."""
    login(client)
    content = (HEADER_ROW + _row("", "Caf\xe9", "Shelf", "1", "g")).encode("latin-1")
    response = _post(client, content, preview=True)
    assert response.status_code == 400
    assert "UTF-8" in response.json["message"]


def test_missing_columns_are_reported(client: FlaskClient, import_setup, sniffed_type):
    """The required columns must be present."""
    login(client)
    response = _post(client, b"Name,Location\nWater,Shelf\n", preview=True)
    assert response.status_code == 400
    assert "Quantity" in response.json["message"]


def test_too_many_rows_are_refused(client: FlaskClient, import_setup, sniffed_type):
    """Files over the row limit are refused."""
    login(client)
    content = (HEADER_ROW + KNOWN_ROW * (MAX_ROWS + 1)).encode()
    response = _post(client, content, preview=True)
    assert response.status_code == 400
    assert "at most" in response.json["message"]


# --- Importing ---------------------------------------------------------------


def test_preview_then_import(client: FlaskClient, import_setup, sniffed_type):
    """
    A preview reports what would happen and saves nothing; the import then
    adds the containers, the new location and the new novel compound. The
    file has Excel's byte order mark and a unit typed in lower case.
    """
    login(client)
    content = "﻿" + (
        HEADER_ROW
        + _row(
            COMPOUND_CAS,
            COMPOUND_NAME,
            "Import cupboard",
            "5",
            "g",
            hazard_codes="H901",
            location_description="Flammables",
        )
        + _row(
            "",
            "Importium",
            "Import cupboard",
            "250",
            "ml",
            smiles=HEXANEDIOL,
            hazard_codes="H900",
        )
        + _row("", "Importium", "Import cupboard", "100", "mL", smiles=HEXANEDIOL)
    )
    workbook_id = import_setup["workbook_id"]

    preview = _post(client, content.encode(), preview=True, workbook_id=workbook_id)
    assert preview.status_code == 200
    assert preview.json["saved"] is False
    assert preview.json["summary"] == {
        "rowCount": 3,
        "errorCount": 0,
        "duplicateCount": 0,
        "newItemCount": 3,
        "newLocationCount": 1,
        "newCompoundCount": 1,
    }
    rows = preview.json["rows"]
    assert [row["row"] for row in rows] == [2, 3, 4]
    assert rows[0]["compoundStatus"] == "existing"
    assert rows[0]["locationStatus"] == "new"
    assert rows[0]["warnings"]  # hazard codes for a known compound are ignored
    assert rows[1]["compoundStatus"] == "new_novel"
    # the second row reuses the compound and location the first one made
    assert rows[2]["compoundStatus"] == "new_novel"
    assert rows[2]["locationStatus"] == "new"
    assert _containers_of(client, "Importium") == []

    response = _post(client, content.encode(), workbook_id=workbook_id)
    assert response.status_code == 201
    assert response.json["saved"] is True

    inventory = _inventory(client)
    location = next(
        loc for loc in inventory["locations"] if loc["name"] == "Import cupboard"
    )
    assert location["description"] == "Flammables"
    novel_items = _containers_of(client, "Importium")
    assert sorted((item["quantity"], item["unit"]) for item in novel_items) == [
        (100.0, "mL"),
        (250.0, "mL"),
    ]
    assert len(_containers_of(client, COMPOUND_NAME)) == 1


def test_import_with_a_bad_row_saves_nothing(
    client: FlaskClient, import_setup, sniffed_type
):
    """One bad row stops the whole import, and its row number is reported."""
    login(client)
    content = (
        HEADER_ROW
        + _row(COMPOUND_CAS, COMPOUND_NAME, "Rollback shelf", "1", "g")
        + _row(COMPOUND_CAS, COMPOUND_NAME, "Rollback shelf", "lots", "g")
        + _row(COMPOUND_CAS, COMPOUND_NAME, "Rollback shelf", "1", "bucket")
    ).encode()
    response = _post(client, content)
    assert response.status_code == 400
    assert response.json["saved"] is False
    errors = {row["row"]: row["error"] for row in response.json["rows"]}
    assert errors[2] is None
    assert "number" in errors[3]
    assert "Unit" in errors[4]
    names = [loc["name"] for loc in _inventory(client)["locations"]]
    assert "Rollback shelf" not in names


def test_unknown_compound_needs_a_workbook(
    client: FlaskClient, import_setup, sniffed_type
):
    """A compound that isn't in the database can't be added without a workbook."""
    login(client)
    content = (
        HEADER_ROW + _row("", "Unknownium", "Shelf", "1", "g", smiles=BUTANEDIOL)
    ).encode()
    response = _post(client, content, preview=True)
    assert response.status_code == 200
    assert "workbook" in response.json["rows"][0]["error"]


def test_new_compound_needs_valid_smiles(
    client: FlaskClient, import_setup, sniffed_type
):
    """New compounds need a SMILES, and it must be one RDKit can read."""
    login(client)
    content = (
        HEADER_ROW
        + _row("", "Unknownium", "Shelf", "1", "g")
        + _row("", "Unknownium", "Shelf", "1", "g", smiles="C1CC(")
    ).encode()
    response = _post(
        client, content, preview=True, workbook_id=import_setup["workbook_id"]
    )
    assert response.status_code == 200
    errors = [row["error"] for row in response.json["rows"]]
    assert "Add its SMILES" in errors[0]
    assert "SMILES isn't valid" in errors[1]


def test_smiles_matches_compound_by_structure(
    client: FlaskClient, import_setup, sniffed_type
):
    """A row whose SMILES is a compound already made under another name uses
    that compound, with a warning, rather than adding a duplicate."""
    login(client)
    content = (
        HEADER_ROW
        + _row("", "Smilesium", "Shelf", "1", "g", smiles=BUTANEDIOL)
        + _row("", "Tetramethylene glycol", "Shelf", "1", "g", smiles=BUTANEDIOL)
    ).encode()
    response = _post(
        client, content, preview=True, workbook_id=import_setup["workbook_id"]
    )
    assert response.status_code == 200
    second = response.json["rows"][1]
    assert second["error"] is None
    assert second["compound"] == "Smilesium"
    assert second["compoundStatus"] == "new_novel"
    assert "SMILES matched Smilesium" in second["warnings"][0]
    assert response.json["summary"]["newCompoundCount"] == 1


def test_workbook_must_be_the_users(client: FlaskClient, import_setup, sniffed_type):
    """Novel compounds can only go into the user's own workbooks."""
    login(client)
    content = (
        HEADER_ROW + _row("", "Unknownium", "Shelf", "1", "g", smiles=BUTANEDIOL)
    ).encode()
    response = _post(client, content, preview=True, workbook_id=999999)
    assert response.status_code == 404


# --- Duplicates --------------------------------------------------------------

AMINE_CAS = "123-45-7"  # Testamine, seeded in tests/populate_database.py
BARCODE_HEADER = HEADER_ROW.rstrip("\n") + ",Barcode\n"


def _barcode_row(cas: str, name: str, location: str, barcode: str) -> str:
    """A row of 1 g of a compound, with a barcode column."""
    return _row(cas, name, location, "1", "g").rstrip("\n") + "," + barcode + "\n"


def _amine_row(location: str, barcode: str = "") -> str:
    """A row of 1 g of Testamine."""
    return _barcode_row(AMINE_CAS, "Testamine", location, barcode)


def test_reimporting_a_file_skips_barcoded_containers(
    client: FlaskClient, import_setup, sniffed_type
):
    """
    Uploading a file again skips the rows whose barcodes are already in the
    inventory, even if the container has moved. Rows without a barcode are
    added again, as identical bottles can share every other detail.
    """
    login(client)
    content = (
        BARCODE_HEADER
        + _amine_row("Dupe shelf", "BC-001")
        + _amine_row("Dupe shelf", "BC-002")
        + _amine_row("Dupe shelf")
    ).encode()
    assert _post(client, content).status_code == 201

    again = _post(client, content)
    assert again.status_code == 201
    summary = again.json["summary"]
    assert (summary["duplicateCount"], summary["newItemCount"]) == (2, 1)
    assert [row["duplicate"] for row in again.json["rows"]] == [True, True, False]
    assert len(_containers_of(client, "Testamine")) == 4

    moved = (BARCODE_HEADER + _amine_row("Moved shelf", "BC-001")).encode()
    preview = _post(client, moved, preview=True)
    assert preview.json["rows"][0]["duplicate"] is True
    # the skipped row doesn't create its location
    assert preview.json["summary"]["newLocationCount"] == 0


def test_barcode_conflicts_are_errors(client: FlaskClient, import_setup, sniffed_type):
    """A barcode on two rows, or already on a different compound's container,
    is reported rather than skipped, as one of them must be wrong."""
    login(client)
    content = (
        BARCODE_HEADER
        + _amine_row("Clash shelf", "BC-100")
        + _amine_row("Clash shelf", "BC-100")
        + _barcode_row(COMPOUND_CAS, COMPOUND_NAME, "Clash shelf", "BC-001")
    ).encode()
    response = _post(client, content, preview=True)
    errors = [row["error"] for row in response.json["rows"]]
    assert errors[0] is None
    assert "earlier row" in errors[1]
    assert "Testamine" in errors[2]
