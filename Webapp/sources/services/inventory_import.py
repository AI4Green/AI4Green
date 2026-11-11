"""
Importing a workgroup's stock from a CSV file, one container per row.

How it works:
- The upload is validated (extension, MIME type sniffed from its bytes, size)
  before anything reads it, then decoded strictly as UTF-8.
- Every row goes through the same service functions as the inventory modals
  (services.inventory.create_item, create_location, create_novel_compound),
  so the import can't accept anything the modals would refuse.
- Compounds are matched by CAS number, then by structure (the SMILES's
  InChI), then by name. Compounds that aren't found become novel compounds,
  which need a SMILES so their formula, molecular weight and InChI can be
  worked out and checked for duplicates.
- Rows whose barcode is already on one of the workgroup's containers are
  skipped, so uploading the same file twice doesn't add everything again.
  Rows without a barcode are always added, as identical bottles can share
  every other detail.
- Each row runs in its own savepoint, so a bad row is reported and undone
  without stopping the rest being checked. Nothing is committed unless every
  row is valid, so an import is saved whole or not at all. A preview runs the
  same steps and always rolls back.

Security notes:
- Nothing in the file is read until validate_upload has passed.
- Row, column and cell counts are capped, as well as the file size.
- Only principal investigators and senior researchers reach this (the API
  checks), and novel compounds can only go into the user's own workbooks.
"""

import csv
import io
import logging
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from flask import current_app
from sources import models, services
from sources.extensions import db
from sources.services.inventory import (
    ALLOWED_UNITS,
    CAS_PATTERN,
    Chemical,
    InventoryNotFoundError,
    InventoryValidationError,
)
from werkzeug.datastructures import FileStorage

logger = logging.getLogger(__name__)

MAX_FILE_BYTES = 1024 * 1024
MAX_ROWS = 2000
MAX_COLUMNS = 30
MAX_CELL_LENGTH = 1000
# CSV types libmagic reports. Only those also in the app's UPLOAD_MIME_TYPES
# are accepted, so the config stays the overall allowlist.
CSV_MIME_TYPES = ("text/csv", "text/plain")
CSV_EXTENSION = ".csv"

FILE_REJECTED_MESSAGE = (
    f"Upload a CSV file (.csv) of up to {MAX_FILE_BYTES // 1024 // 1024} MB."
)
ENCODING_MESSAGE = "Save the file as CSV UTF-8 and upload it again."

# Column names users may use, after _normalise_header, for each field.
COLUMN_ALIASES: Dict[str, tuple] = {
    "cas": ("cas", "cas number", "cas no"),
    "name": ("name", "compound", "compound name"),
    "smiles": ("smiles", "smiles string"),
    "hazard_codes": ("hazard codes", "hazards", "h codes"),
    "location": ("location", "location name"),
    "location_description": (
        "location description",
        "location summary",
    ),
    "quantity": ("quantity", "amount"),
    "unit": ("unit", "units"),
    "supplier": ("supplier",),
    "lot_number": ("lot number", "lot"),
    "purity": ("purity",),
    "barcode": ("barcode",),
    "date_received": ("date received",),
    "expiry_date": ("expiry date",),
}
REQUIRED_COLUMNS = ("name", "location", "quantity", "unit")
# Units as users often type them, mapped to the allowed units.
UNIT_SPELLINGS = {
    **{unit.lower(): unit for unit in ALLOWED_UNITS},
    "ul": "µL",
    "μl": "µL",  # Greek mu, which looks the same as the micro sign
}


@dataclass
class RowResult:
    """What importing one row did, or why it couldn't."""

    row: int
    name: str
    compound: Optional[str] = None
    compound_status: Optional[str] = None  # existing | new_novel
    location: Optional[str] = None
    location_status: Optional[str] = None  # existing | new
    quantity: Optional[str] = None
    unit: Optional[str] = None
    # the container is already in the inventory, so the row is skipped
    duplicate: bool = False
    error: Optional[str] = None
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """The result as the front end uses it."""
        return {
            "row": self.row,
            "name": self.name,
            "compound": self.compound,
            "compoundStatus": self.compound_status,
            "location": self.location,
            "locationStatus": self.location_status,
            "quantity": self.quantity,
            "unit": self.unit,
            "duplicate": self.duplicate,
            "error": self.error,
            "warnings": self.warnings,
        }


# --- Reading the file --------------------------------------------------------


def _allowed_mime_types() -> List[str]:
    """The CSV MIME types that the app-wide upload allowlist also permits."""
    configured = current_app.config["UPLOAD_MIME_TYPES"]
    return [mime_type for mime_type in CSV_MIME_TYPES if mime_type in configured]


def _read_validated_text(file: Optional[FileStorage]) -> str:
    """
    Validates the upload and only then reads and decodes it.

    Args:
        file: The uploaded file, if there is one.

    Returns:
        The file's text.

    Raises:
        InventoryValidationError: If there's no file, it fails validation, or
            it isn't UTF-8.
    """
    if file is None:
        raise InventoryValidationError("Choose a CSV file to import.")
    if CSV_EXTENSION not in current_app.config["UPLOAD_EXTENSIONS"]:
        raise InventoryValidationError(FILE_REJECTED_MESSAGE)
    error = services.file_attachments.validate_upload(
        file,
        allowed_extensions=[CSV_EXTENSION],
        allowed_mime_types=_allowed_mime_types(),
        max_bytes=MAX_FILE_BYTES,
    )
    if error is not None:
        raise InventoryValidationError(FILE_REJECTED_MESSAGE)

    content = file.stream.read(MAX_FILE_BYTES + 1)
    try:
        # utf-8-sig drops the byte order mark Excel adds to "CSV UTF-8" files
        return content.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise InventoryValidationError(ENCODING_MESSAGE) from None


def _normalise_header(header: str) -> str:
    """Lower case, with runs of anything but letters and digits as one space."""
    return re.sub(r"[^a-z0-9]+", " ", header.lower()).strip()


def _map_columns(headers: List[str]) -> Dict[str, int]:
    """
    Works out which column holds each field.

    Args:
        headers: The header row.

    Returns:
        Each field's column index.

    Raises:
        InventoryValidationError: If a field appears twice or a required one
            is missing.
    """
    field_by_alias = {
        alias: field_name
        for field_name, aliases in COLUMN_ALIASES.items()
        for alias in aliases
    }
    columns: Dict[str, int] = {}
    for index, header in enumerate(headers):
        field_name = field_by_alias.get(_normalise_header(header))
        if field_name is None:
            continue
        if field_name in columns:
            raise InventoryValidationError(
                f"The file has more than one {field_name.replace('_', ' ')} column."
            )
        columns[field_name] = index

    missing = [name for name in REQUIRED_COLUMNS if name not in columns]
    if missing:
        raise InventoryValidationError(
            "The file needs these columns: "
            + ", ".join(name.replace("_", " ").capitalize() for name in missing)
            + "."
        )
    return columns


def _read_rows(text: str) -> List[Dict[str, Any]]:
    """
    Parses the CSV text into rows keyed by field, skipping blank rows.

    Args:
        text: The decoded file.

    Returns:
        Rows as `{"row": spreadsheet row number, "cells": {field: text}}`.

    Raises:
        InventoryValidationError: If the file is empty, malformed, or over
            the row, column or cell limits.
    """
    reader = csv.reader(io.StringIO(text, newline=""))
    try:
        headers = next(reader, None)
        if headers is None:
            raise InventoryValidationError("The file is empty.")
        if len(headers) > MAX_COLUMNS:
            raise InventoryValidationError(
                f"The file can have at most {MAX_COLUMNS} columns."
            )
        columns = _map_columns(headers)

        rows = []
        for values in reader:
            if not any(value.strip() for value in values):
                continue
            # the header is spreadsheet row 1
            row_number = reader.line_num
            if len(rows) >= MAX_ROWS:
                raise InventoryValidationError(
                    f"The file can have at most {MAX_ROWS} rows. Split it into "
                    "smaller files."
                )
            if len(values) > MAX_COLUMNS:
                raise InventoryValidationError(
                    f"Row {row_number} has more than {MAX_COLUMNS} columns."
                )
            if any(len(value) > MAX_CELL_LENGTH for value in values):
                raise InventoryValidationError(
                    f"Row {row_number} has a value over {MAX_CELL_LENGTH} "
                    "characters."
                )
            cells = {
                field_name: values[index].strip() if index < len(values) else ""
                for field_name, index in columns.items()
            }
            rows.append({"row": row_number, "cells": cells})
    except csv.Error:
        raise InventoryValidationError(
            "The file couldn't be read as CSV. Check it was saved as CSV."
        ) from None

    if not rows:
        raise InventoryValidationError("The file has no rows to import.")
    return rows


# --- Importing rows ----------------------------------------------------------


@dataclass
class _ImportContext:
    """
    State shared by the rows of one import. The id sets only hold locations
    and compounds from rows that succeeded, as a failed row's are rolled back.
    """

    workgroup: models.WorkGroup
    workbook: Optional[models.WorkBook]
    new_location_ids: set = field(default_factory=set)
    new_compound_ids: set = field(default_factory=set)
    # barcodes on earlier rows, so the file can't give one to two containers
    barcodes: set = field(default_factory=set)
    new_item_count: int = 0
    # made by the row being imported, kept only if the row succeeds
    pending_location_id: Optional[int] = None
    pending_compound_id: Optional[int] = None

    def keep_pending(self) -> None:
        """Records the current row's new location and compound as kept."""
        if self.pending_location_id is not None:
            self.new_location_ids.add(self.pending_location_id)
        if self.pending_compound_id is not None:
            self.new_compound_ids.add(self.pending_compound_id)
        self.clear_pending()

    def clear_pending(self) -> None:
        """Forgets the current row's new location and compound."""
        self.pending_location_id = None
        self.pending_compound_id = None


def _read_quantity_text(text: str) -> float:
    """
    Reads a quantity cell as a number, which create_item then validates.

    Args:
        text: The cell.

    Returns:
        The number.
    """
    if not text:
        raise InventoryValidationError("Enter a quantity.")
    try:
        return float(text)
    except ValueError:
        raise InventoryValidationError("Quantity must be a number.") from None


def _resolve_location(
    context: _ImportContext, cells: Dict[str, str], result: RowResult
) -> models.Location:
    """
    Finds the row's location by name, or creates it with the row's
    description.

    Args:
        context: The import.
        cells: The row.
        result: The row's result, updated with the location.

    Returns:
        The location.
    """
    name = cells["location"]
    if not name:
        raise InventoryValidationError("Enter a location.")
    location = services.inventory.find_location(context.workgroup, name)
    if location is None:
        location = services.inventory.create_location(
            context.workgroup,
            {"name": name, "description": cells.get("location_description")},
        )
        context.pending_location_id = location.id
    result.location = location.name
    is_new = location.id in (context.new_location_ids | {context.pending_location_id})
    result.location_status = "new" if is_new else "existing"
    return location


def _read_inchi(smiles: str) -> Optional[str]:
    """
    Works out the structure's InChI from the row's SMILES, if it has one.

    Args:
        smiles: The SMILES cell, which may be blank.

    Returns:
        The InChI, or None if the cell is blank.

    Raises:
        InventoryValidationError: If RDKit can't read the SMILES.
    """
    if not smiles:
        return None
    inchi = services.all_compounds.smiles_to_inchi(smiles)
    if not inchi:
        raise InventoryValidationError(
            "The SMILES isn't valid. Check it, e.g. by drawing it in the sketcher."
        )
    return inchi


def _find_chemical(
    context: _ImportContext,
    cas: str,
    inchi: Optional[str],
    name: str,
    result: RowResult,
) -> Optional[Chemical]:
    """
    Finds the row's compound by CAS number, then by structure, then (only for
    rows without a CAS number) by name, warning about anything that doesn't
    agree so the user can check it in the preview.

    Args:
        context: The import.
        cas: The CAS number, or blank.
        inchi: The InChI of the row's SMILES, or None.
        name: The compound name.
        result: The row's result, given any warnings.

    Returns:
        The compound or novel compound, or None if it isn't in the database.
    """
    workgroup = context.workgroup
    if cas:
        chemical = services.inventory.find_chemical(workgroup, cas)
        if chemical is not None:
            if chemical.name.lower() != name.lower():
                # often a synonym, but could be a typo in the CAS number
                result.warnings.append(f"CAS {cas} matched {chemical.name}.")
            if inchi and chemical.inchi and chemical.inchi != inchi:
                result.warnings.append(
                    f"The SMILES isn't the structure of {chemical.name}, which "
                    f"CAS {cas} matched."
                )
            return chemical

    if inchi:
        chemical = services.inventory.find_chemical_by_inchi(workgroup, inchi)
        if chemical is not None:
            if cas:
                result.warnings.append(
                    f"CAS {cas} isn't in the database, but the SMILES matched Novel Compound "
                    f"{chemical.name}."
                )
            elif chemical.name.lower() != name.lower():
                result.warnings.append(f"The SMILES matched {chemical.name}.")
            return chemical

    if cas:
        return None
    return services.inventory.find_chemical(workgroup, name)


def _resolve_chemical(
    context: _ImportContext, cells: Dict[str, str], result: RowResult
) -> Chemical:
    """
    Finds the row's compound, or adds it as a novel compound to the chosen
    workbook. New compounds need a SMILES, from which their formula,
    molecular weight and InChI are worked out.

    Args:
        context: The import.
        cells: The row.
        result: The row's result, updated with the compound.

    Returns:
        The compound or novel compound.
    """
    cas = cells.get("cas", "")
    name = cells["name"]
    smiles = cells.get("smiles", "")
    hazard_codes = cells.get("hazard_codes", "")
    if not name:
        raise InventoryValidationError("Enter a name.")
    if cas and not CAS_PATTERN.fullmatch(cas):
        raise InventoryValidationError(
            f"{cas} isn't a CAS number. Use the form 67-64-1."
        )
    inchi = _read_inchi(smiles)

    chemical = _find_chemical(context, cas, inchi, name, result)
    if chemical is None:
        if context.workbook is None:
            raise InventoryValidationError(
                "This compound isn't in the database. Choose a workbook to "
                "add new compounds to."
            )
        if not smiles:
            raise InventoryValidationError(
                "This compound isn't in the database. Add its SMILES so it "
                "can be added as a new compound."
            )
        chemical = services.inventory.create_novel_compound(
            context.workbook,
            {
                "name": name,
                "cas": cas or None,
                "smiles": smiles,
                "hazardCodes": hazard_codes or None,
            },
        )
        context.pending_compound_id = chemical.id
    elif hazard_codes and isinstance(chemical, models.Compound):
        result.warnings.append(
            "Hazard codes are taken from the compound database, not the file."
        )

    is_new = isinstance(chemical, models.NovelCompound) and chemical.id in (
        context.new_compound_ids | {context.pending_compound_id}
    )
    result.compound = chemical.name
    result.compound_status = "new_novel" if is_new else "existing"
    return chemical


def _is_already_held(context: _ImportContext, barcode: str, chemical: Chemical) -> bool:
    """
    Checks whether the row's container is already in the inventory, which is
    judged by barcode alone: rows without one are always added.

    Args:
        context: The import.
        barcode: The row's barcode, or blank.
        chemical: The row's compound.

    Returns:
        True if a container of this compound already has the barcode. If a
        different compound's container has it, this returns False and
        create_item refuses the barcode, as one of them must be wrong.

    Raises:
        InventoryValidationError: If an earlier row has the barcode too.
    """
    if not barcode:
        return False
    if barcode in context.barcodes:
        raise InventoryValidationError(f"Barcode {barcode} is on an earlier row.")
    context.barcodes.add(barcode)

    existing = services.inventory.find_item_by_barcode(context.workgroup, barcode)
    if existing is None:
        return False
    held = existing.workgroup_compound.chemical
    return type(held) is type(chemical) and held.id == chemical.id


def _import_row(
    context: _ImportContext, row_number: int, cells: Dict[str, str]
) -> RowResult:
    """
    Imports one row in a savepoint, which is rolled back if the row fails.

    Args:
        context: The import.
        row_number: The spreadsheet row number, for messages.
        cells: The row.

    Returns:
        What happened, including the error if it failed.
    """
    result = RowResult(
        row=row_number,
        name=cells["name"],
        quantity=cells["quantity"],
        unit=cells["unit"],
    )
    savepoint = db.session.begin_nested()
    try:
        location = _resolve_location(context, cells, result)
        chemical = _resolve_chemical(context, cells, result)
        unit_text = cells["unit"]
        item_data = {
            "compound": {
                "source": (
                    "novel"
                    if isinstance(chemical, models.NovelCompound)
                    else "compound"
                ),
                "id": chemical.id,
            },
            "locationId": location.id,
            "quantity": _read_quantity_text(cells["quantity"]),
            "unit": UNIT_SPELLINGS.get(unit_text.lower(), unit_text),
            "supplier": cells.get("supplier"),
            "lotNumber": cells.get("lot_number"),
            "purity": cells.get("purity"),
            "barcode": cells.get("barcode"),
            "dateReceived": cells.get("date_received"),
            "expiryDate": cells.get("expiry_date"),
        }
        if _is_already_held(context, cells.get("barcode", ""), chemical):
            # undo anything the row created, e.g. a new location
            savepoint.rollback()
            context.clear_pending()
            result.duplicate = True
            result.warnings.append(
                "Already in the inventory, so this row will be skipped."
            )
            return result

        services.inventory.create_item(context.workgroup, item_data)
        savepoint.commit()
        context.keep_pending()
        context.new_item_count += 1
    except InventoryValidationError as error:
        savepoint.rollback()
        context.clear_pending()
        result.error = str(error)
    except InventoryNotFoundError:
        savepoint.rollback()
        context.clear_pending()
        result.error = "The compound or location couldn't be found."
    except ValueError:
        # e.g. the ledger refusing a change; its message isn't for users
        savepoint.rollback()
        context.clear_pending()
        logger.exception("Inventory import row %s failed", row_number)
        result.error = "This row couldn't be imported."
    return result


def import_csv(
    workgroup: models.WorkGroup,
    file: Optional[FileStorage],
    workbook_id: Any,
    save: bool,
) -> Dict[str, Any]:
    """
    Validates a CSV upload and imports its rows, one container per row,
    creating any locations and novel compounds they need.

    Args:
        workgroup: The workgroup.
        file: The uploaded CSV file.
        workbook_id: The id of the user's workbook for new novel compounds,
            as text from the form, or None/blank if there isn't one.
        save: True to save the import if every row is valid. False to only
            preview it; nothing is saved either way.

    Returns:
        `saved`, the per-row `rows` (see RowResult.to_dict) and a `summary`
        of the counts.

    Raises:
        InventoryValidationError: If the file, its columns or the workbook
            id are invalid. Row problems are reported in `rows` instead.
        InventoryNotFoundError: If the workbook isn't one of the user's.
    """
    text = _read_validated_text(file)
    rows = _read_rows(text)
    workbook = None
    if workbook_id not in (None, ""):
        try:
            workbook_id = int(workbook_id)
        except (TypeError, ValueError):
            raise InventoryValidationError("Choose a valid workbook.") from None
        workbook = services.inventory.get_user_workbook(workgroup, workbook_id)

    context = _ImportContext(workgroup=workgroup, workbook=workbook)
    try:
        results = [_import_row(context, row["row"], row["cells"]) for row in rows]
        error_count = sum(1 for result in results if result.error)
        duplicate_count = sum(1 for result in results if result.duplicate)
        saved = save and error_count == 0
        if saved:
            db.session.commit()
            logger.info(
                "Imported %s containers into workgroup %s, skipping %s duplicates",
                context.new_item_count,
                workgroup.id,
                duplicate_count,
            )
        else:
            db.session.rollback()
    except Exception:
        db.session.rollback()
        raise

    return {
        "saved": saved,
        "rows": [result.to_dict() for result in results],
        "summary": {
            "rowCount": len(results),
            "errorCount": error_count,
            "duplicateCount": duplicate_count,
            # containers the import adds (or would add, for a preview)
            "newItemCount": context.new_item_count,
            "newLocationCount": len(context.new_location_ids),
            "newCompoundCount": len(context.new_compound_ids),
        },
    }
