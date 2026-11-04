"""
Inventory services: reading a workgroup's stock and making changes to it.

Security notes:
- Every lookup by id also filters by the workgroup, so an id from another
  workgroup is treated as not found and can't be read or changed.
- Novel compounds are only visible through the workbooks the current user is
  a member of, matching how the reaction table treats them.
- Input is validated here, whatever the front end has already checked.
  Problems raise InventoryValidationError with a message that is safe to show
  to the user.
"""

import math
import re
from collections import defaultdict
from datetime import date
from decimal import Decimal
from typing import Any, Dict, List, Optional, Union

from flask_login import current_user
from sources import models, services
from sources.extensions import db
from sources.models.inventory_ledger import StockChangeReason
from sqlalchemy import func
from sqlalchemy.orm import joinedload

# Roles allowed to change the inventory. Any workgroup member can read it.
EDITOR_ROLES = ("principal_investigator", "senior_researcher")

EXPIRING_SOON_DAYS = 30
ALLOWED_UNITS = ("mg", "g", "kg", "µL", "mL", "L")
# the largest value a Numeric(12, 4) column holds
MAX_QUANTITY = Decimal("99999999.9999")
QUANTITY_PLACES = Decimal("0.0001")
NOTES_MAX_LENGTH = 500
LOCATION_NAME_MAX_LENGTH = 100
LOCATION_DESCRIPTION_MAX_LENGTH = 255
MAX_LOOKUP_QUERY_LENGTH = 200
MAX_SUGGESTIONS = 5

# digits-digits-check digit, without a leading zero, as the reaction table uses
CAS_PATTERN = re.compile(r"[1-9]\d{1,6}-\d{2}-\d")
HAZARD_CODE_PATTERN = re.compile(r"(EU)?H\d{3}[A-Za-z]{0,2}")

Chemical = Union[models.Compound, models.NovelCompound]


class InventoryValidationError(ValueError):
    """Invalid input. The message is shown to the user."""


class InventoryNotFoundError(LookupError):
    """A requested row doesn't exist in this workgroup, or the user can't see it."""


# --- Access ------------------------------------------------------------------


def can_edit(workgroup: models.WorkGroup) -> bool:
    """
    Whether the current user may change the workgroup's inventory.

    Args:
        workgroup: The workgroup.

    Returns:
        True for principal investigators and senior researchers.
    """
    role = services.workgroup.get_user_type(workgroup.name, current_user)
    return role in EDITOR_ROLES


# --- Input validation --------------------------------------------------------


def _read_id(value: Any, label: str) -> int:
    """Reads a required positive integer id."""
    # bool is a subclass of int, so rule it out explicitly
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise InventoryValidationError(f"Choose a valid {label.lower()}.")
    return value


def _read_text(value: Any, label: str, max_length: int) -> Optional[str]:
    """Reads optional text: trimmed, with blank or missing values as None."""
    if value is None:
        return None
    if not isinstance(value, str):
        raise InventoryValidationError(f"{label} must be text.")
    text = value.strip()
    if len(text) > max_length:
        raise InventoryValidationError(
            f"{label} must be {max_length} characters or fewer."
        )
    return text or None


def _read_required_text(value: Any, label: str, max_length: int) -> str:
    """Reads required text, trimmed."""
    text = _read_text(value, label, max_length)
    if text is None:
        raise InventoryValidationError(f"Enter a {label.lower()}.")
    return text


def _read_number(value: Any, label: str) -> float:
    """Reads a finite JSON number, rejecting booleans, strings, NaN and Infinity."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InventoryValidationError(f"{label} must be a number.")
    if not math.isfinite(value):
        raise InventoryValidationError(f"{label} must be a finite number.")
    return value


def _read_quantity(value: Any, label: str, allow_zero: bool) -> Decimal:
    """Reads a quantity, rounded to the 4 decimal places the database stores."""
    number = _read_number(value, label)
    quantity = Decimal(str(number)).quantize(QUANTITY_PLACES)
    if quantity < 0 or (quantity == 0 and not allow_zero):
        smallest = "zero or more" if allow_zero else "more than zero"
        raise InventoryValidationError(f"{label} must be {smallest}.")
    if quantity > MAX_QUANTITY:
        raise InventoryValidationError(f"{label} is too large.")
    return quantity


def _read_optional_positive(value: Any, label: str, maximum: float) -> Optional[float]:
    """Reads an optional positive number up to `maximum`; blank is None."""
    if value is None or value == "":
        return None
    number = _read_number(value, label)
    if number <= 0 or number > maximum:
        raise InventoryValidationError(
            f"{label} must be more than zero and at most {maximum:g}."  # noqa: E231
        )
    return float(number)


def _read_date(value: Any, label: str) -> Optional[date]:
    """Reads an optional YYYY-MM-DD date; blank is None."""
    if value is None or value == "":
        return None
    if not isinstance(value, str):
        raise InventoryValidationError(f"{label} must be a date.")
    try:
        return date.fromisoformat(value)
    except ValueError:
        raise InventoryValidationError(f"{label} must be a valid date.") from None


def _read_unit(value: Any) -> str:
    """Reads a unit from the allowed list."""
    if value not in ALLOWED_UNITS:
        raise InventoryValidationError(
            f"Unit must be one of {', '.join(ALLOWED_UNITS)}."
        )
    return value


# --- Lookups scoped to a workgroup ------------------------------------------


def get_item(workgroup: models.WorkGroup, item_id: int) -> models.InventoryItem:
    """
    Gets one of the workgroup's containers.

    Args:
        workgroup: The workgroup.
        item_id: The InventoryItem id.

    Returns:
        The container.

    Raises:
        InventoryNotFoundError: If it doesn't exist in this workgroup.
    """
    item = (
        db.session.query(models.InventoryItem)
        .join(models.WorkgroupCompound)
        .filter(
            models.InventoryItem.id == item_id,
            models.WorkgroupCompound.workgroup_id == workgroup.id,
        )
        .first()
    )
    if item is None:
        raise InventoryNotFoundError()
    return item


def _get_active_item(workgroup: models.WorkGroup, item_id: int) -> models.InventoryItem:
    """Gets a container that hasn't been archived (disposed of)."""
    item = get_item(workgroup, item_id)
    if item.is_archived:
        raise InventoryValidationError("This container has already been disposed of.")
    return item


def _get_location(workgroup: models.WorkGroup, location_id: int) -> models.Location:
    """Gets one of the workgroup's locations, or raises InventoryNotFoundError."""
    location = (
        db.session.query(models.Location)
        .filter(
            models.Location.id == location_id,
            models.Location.workgroup_id == workgroup.id,
        )
        .first()
    )
    if location is None:
        raise InventoryNotFoundError()
    return location


def _get_active_location(
    workgroup: models.WorkGroup, location_id: int
) -> models.Location:
    """Gets a location that containers can be put in."""
    location = _get_location(workgroup, location_id)
    if location.is_archived:
        raise InventoryValidationError(
            "This location is archived. Restore it before storing containers there."
        )
    return location


def _get_user_novel_compound(
    workgroup: models.WorkGroup, novel_compound_id: int
) -> models.NovelCompound:
    """
    Gets a novel compound from one of the user's workbooks in this workgroup,
    or raises InventoryNotFoundError.
    """
    workbook_ids = [
        wb.id
        for wb in services.workbook.get_workbooks_from_user_group_combination(
            workgroup.name
        )
    ]
    novel_compound = (
        db.session.query(models.NovelCompound)
        .filter(
            models.NovelCompound.id == novel_compound_id,
            models.NovelCompound.workbook.in_(workbook_ids),
        )
        .first()
    )
    if novel_compound is None:
        raise InventoryNotFoundError()
    return novel_compound


def _find_workgroup_compound(
    workgroup: models.WorkGroup, chemical: Chemical
) -> Optional[models.WorkgroupCompound]:
    """Finds the workgroup's row for a compound or novel compound, if it has one."""
    query = db.session.query(models.WorkgroupCompound).filter(
        models.WorkgroupCompound.workgroup_id == workgroup.id
    )
    if isinstance(chemical, models.NovelCompound):
        query = query.filter(models.WorkgroupCompound.novel_compound_id == chemical.id)
    else:
        query = query.filter(models.WorkgroupCompound.compound_id == chemical.id)
    return query.first()


# --- Serialisation -----------------------------------------------------------


def _hazard_codes(hphrase: Optional[str]) -> List[str]:
    """Splits a hazard phrase into its codes, dropping placeholders such as
    "Unknown"."""
    codes = (
        code.strip() for code in services.hazard_code.string_to_list(hphrase or "")
    )
    return [code for code in codes if HAZARD_CODE_PATTERN.fullmatch(code)]


def _iso(value: Optional[date]) -> Optional[str]:
    """Formats an optional date or datetime as ISO 8601."""
    return value.isoformat() if value else None


def _number(value: Optional[Decimal]) -> Optional[float]:
    """Converts an optional Decimal to a float for JSON."""
    return float(value) if value is not None else None


def serialise_location(location: models.Location) -> Dict[str, Any]:
    """A location as the front end uses it."""
    return {
        "id": location.id,
        "name": location.name,
        "description": location.description,
        "isArchived": location.is_archived,
    }


def serialise_item(item: models.InventoryItem, today: date) -> Dict[str, Any]:
    """A container as the front end uses it."""
    days_to_expiry = (item.expiry_date - today).days if item.expiry_date else None
    return {
        "id": item.id,
        "location": {"id": item.location.id, "name": item.location.name},
        "quantity": _number(item.quantity),
        "unit": item.unit,
        "supplier": item.supplier,
        "lotNumber": item.lot_number,
        "purity": item.purity,
        "barcode": item.barcode,
        "dateReceived": _iso(item.date_received),
        "dateOpened": _iso(item.date_opened),
        "expiryDate": _iso(item.expiry_date),
        "isExpired": days_to_expiry is not None and days_to_expiry < 0,
        "isExpiringSoon": days_to_expiry is not None
        and 0 <= days_to_expiry <= EXPIRING_SOON_DAYS,
        "isArchived": item.is_archived,
    }


def _serialise_workgroup_compound(
    wc: models.WorkgroupCompound, items: List[models.InventoryItem], today: date
) -> Dict[str, Any]:
    """
    A compound the workgroup holds, with its containers. The total only counts
    active containers in the threshold unit, as there's no unit conversion.
    """
    chemical = wc.chemical
    total = sum(
        (
            item.quantity
            for item in items
            if not item.is_archived and item.unit == wc.threshold_unit
        ),
        Decimal(0),
    )
    return {
        "id": wc.id,
        "compound": {
            "name": chemical.name,
            "cas": chemical.cas or "",
            "hazards": _hazard_codes(chemical.hphrase),
            "isNovel": wc.novel_compound_id is not None,
        },
        "preferredSupplier": wc.preferred_supplier,
        "minStockThreshold": _number(wc.min_stock_threshold),
        "thresholdUnit": wc.threshold_unit,
        "totalQuantity": float(total),
        "isBelowThreshold": wc.min_stock_threshold is not None
        and total < wc.min_stock_threshold,
        "notes": wc.notes,
        "items": [serialise_item(item, today) for item in items],
    }


def serialise_chemical(
    chemical: Chemical, workgroup: models.WorkGroup
) -> Dict[str, Any]:
    """
    A compound or novel compound as returned by the lookup, including whether
    the workgroup already holds it.
    """
    is_novel = isinstance(chemical, models.NovelCompound)
    held = _find_workgroup_compound(workgroup, chemical)
    workbook = db.session.get(models.WorkBook, chemical.workbook) if is_novel else None
    return {
        "source": "novel" if is_novel else "compound",
        "id": chemical.id,
        "name": chemical.name,
        "cas": chemical.cas or "",
        "hazards": _hazard_codes(chemical.hphrase),
        "workbook": workbook.name if workbook else None,
        "workgroupCompoundId": held.id if held else None,
    }


def serialise_ledger_entry(entry: models.InventoryLedgerEntry) -> Dict[str, Any]:
    """A ledger entry for the container history."""
    return {
        "id": entry.id,
        "change": float(entry.change),
        "reason": entry.reason.value,
        "notes": entry.notes,
        "createdAt": _iso(entry.created_at),
        "user": {"name": entry.user.fullname} if entry.user else None,
    }


# --- Reading -----------------------------------------------------------------


def get_inventory(workgroup: models.WorkGroup) -> Dict[str, Any]:
    """
    Everything the inventory page needs for a workgroup.

    Args:
        workgroup: The workgroup.

    Returns:
        `workgroup`, `permissions`, `locations`, `workbooks` (the user's, for
        adding novel compounds) and `compounds`, each with its containers.
    """
    today = date.today()

    locations = (
        db.session.query(models.Location)
        .filter(models.Location.workgroup_id == workgroup.id)
        .order_by(models.Location.name)
        .all()
    )
    workgroup_compounds = (
        db.session.query(models.WorkgroupCompound)
        .options(
            joinedload(models.WorkgroupCompound.compound),
            joinedload(models.WorkgroupCompound.novel_compound),
        )
        .filter(models.WorkgroupCompound.workgroup_id == workgroup.id)
        .all()
    )
    items = (
        db.session.query(models.InventoryItem)
        .join(models.WorkgroupCompound)
        .options(joinedload(models.InventoryItem.location))
        .filter(models.WorkgroupCompound.workgroup_id == workgroup.id)
        .order_by(models.InventoryItem.id)
        .all()
    )
    items_by_compound = defaultdict(list)
    for item in items:
        items_by_compound[item.workgroup_compound_id].append(item)

    compounds = [
        _serialise_workgroup_compound(wc, items_by_compound[wc.id], today)
        for wc in workgroup_compounds
    ]
    compounds.sort(key=lambda c: c["compound"]["name"].lower())
    workbooks = sorted(
        services.workbook.get_workbooks_from_user_group_combination(workgroup.name),
        key=lambda wb: wb.name,
    )

    return {
        "workgroup": workgroup.name,
        "permissions": {"canEdit": can_edit(workgroup)},
        "locations": [serialise_location(location) for location in locations],
        "workbooks": [{"id": wb.id, "name": wb.name} for wb in workbooks],
        "compounds": compounds,
    }


def get_item_history(workgroup: models.WorkGroup, item_id: int) -> List[Dict[str, Any]]:
    """
    A container's ledger, oldest first.

    Args:
        workgroup: The workgroup.
        item_id: The InventoryItem id.

    Returns:
        The serialised ledger entries.
    """
    item = get_item(workgroup, item_id)
    entries = item.ledger_entries.order_by(
        models.InventoryLedgerEntry.created_at, models.InventoryLedgerEntry.id
    ).all()
    return [serialise_ledger_entry(entry) for entry in entries]


# --- Changing stock ----------------------------------------------------------


def _resolve_workgroup_compound(
    workgroup: models.WorkGroup, data: Dict[str, Any]
) -> models.WorkgroupCompound:
    """
    Finds the WorkgroupCompound a new container belongs to, creating it the
    first time the workgroup stocks a compound.

    `data` holds either `workgroupCompoundId`, or `compound: {source, id}` from
    the lookup, where source is "compound" or "novel".
    """
    workgroup_compound_id = data.get("workgroupCompoundId")
    compound = data.get("compound")
    if (workgroup_compound_id is None) == (compound is None):
        raise InventoryValidationError("Choose a compound.")

    if workgroup_compound_id is not None:
        wc = (
            db.session.query(models.WorkgroupCompound)
            .filter(
                models.WorkgroupCompound.id
                == _read_id(workgroup_compound_id, "Compound"),
                models.WorkgroupCompound.workgroup_id == workgroup.id,
            )
            .first()
        )
        if wc is None:
            raise InventoryNotFoundError()
        return wc

    if not isinstance(compound, dict) or compound.get("source") not in (
        "compound",
        "novel",
    ):
        raise InventoryValidationError("Choose a compound.")
    compound_id = _read_id(compound.get("id"), "Compound")

    if compound["source"] == "novel":
        chemical = _get_user_novel_compound(workgroup, compound_id)
    else:
        chemical = services.compound.get(compound_id)
        if chemical is None:
            raise InventoryNotFoundError()

    wc = _find_workgroup_compound(workgroup, chemical)
    if wc is None:
        wc = models.WorkgroupCompound(workgroup_id=workgroup.id)
        if isinstance(chemical, models.NovelCompound):
            wc.novel_compound_id = chemical.id
        else:
            wc.compound_id = chemical.id
        db.session.add(wc)
    return wc


def add_item(workgroup: models.WorkGroup, data: Dict[str, Any]) -> models.InventoryItem:
    """
    Adds a container and records it as received in the ledger.

    Args:
        workgroup: The workgroup receiving the stock.
        data: The JSON body from the front end.

    Returns:
        The new container.
    """
    # validate the plain fields before anything is added to the session
    quantity = _read_quantity(data.get("quantity"), "Quantity", allow_zero=False)
    unit = _read_unit(data.get("unit"))
    supplier = _read_text(data.get("supplier"), "Supplier", 255)
    lot_number = _read_text(data.get("lotNumber"), "Lot number", 100)
    purity = _read_text(data.get("purity"), "Purity", 50)
    barcode = _read_text(data.get("barcode"), "Barcode", 100)
    date_received = _read_date(data.get("dateReceived"), "Date received")
    expiry_date = _read_date(data.get("expiryDate"), "Expiry date")
    location = _get_active_location(
        workgroup, _read_id(data.get("locationId"), "Location")
    )

    wc = _resolve_workgroup_compound(workgroup, data)
    # the first container sets the compound's defaults
    if wc.threshold_unit is None:
        wc.threshold_unit = unit
    if wc.preferred_supplier is None and supplier:
        wc.preferred_supplier = supplier

    db.session.flush()  # gives a new WorkgroupCompound its id
    item = services.inventory_ledger.create_inventory_item(
        wc.id,
        location.id,
        quantity,
        unit,
        current_user.id,
        supplier=supplier,
        lot_number=lot_number,
        purity=purity,
        barcode=barcode,
        date_received=date_received,
        expiry_date=expiry_date,
    )
    db.session.commit()
    return item


def adjust_item(
    workgroup: models.WorkGroup, item_id: int, data: Dict[str, Any]
) -> models.InventoryItem:
    """
    Sets a container's quantity to a measured value and records the difference.
    A container whose quantity goes down is marked as opened, if it wasn't.

    Args:
        workgroup: The workgroup.
        item_id: The InventoryItem id.
        data: The JSON body, with `quantity` and optional `notes`.

    Returns:
        The updated container.
    """
    item = _get_active_item(workgroup, item_id)
    quantity = _read_quantity(data.get("quantity"), "Quantity", allow_zero=True)
    notes = _read_text(data.get("notes"), "Notes", NOTES_MAX_LENGTH)

    change = quantity - item.quantity
    if change == 0:
        raise InventoryValidationError("That's the same as the current quantity.")
    if change < 0 and item.date_opened is None:
        item.date_opened = date.today()

    services.inventory_ledger.adjust_stock(
        item.id, change, StockChangeReason.ADJUSTMENT, current_user.id, notes=notes
    )
    db.session.commit()
    return item


def move_item(
    workgroup: models.WorkGroup, item_id: int, data: Dict[str, Any]
) -> models.InventoryItem:
    """
    Moves a container to another of the workgroup's locations. Moves don't
    change the quantity, so they aren't recorded in the ledger.

    Args:
        workgroup: The workgroup.
        item_id: The InventoryItem id.
        data: The JSON body, with `locationId`.

    Returns:
        The updated container.
    """
    item = _get_active_item(workgroup, item_id)
    location = _get_active_location(
        workgroup, _read_id(data.get("locationId"), "Location")
    )
    if location.id == item.location_id:
        raise InventoryValidationError("The container is already in this location.")

    item.location = location
    db.session.commit()
    return item


def dispose_item(
    workgroup: models.WorkGroup, item_id: int, data: Dict[str, Any]
) -> models.InventoryItem:
    """
    Disposes of a container: records the remaining quantity as disposed, sets
    it to zero and archives the container. Its history is kept.

    Args:
        workgroup: The workgroup.
        item_id: The InventoryItem id.
        data: The JSON body, with optional `notes`.

    Returns:
        The archived container.
    """
    item = _get_active_item(workgroup, item_id)
    notes = _read_text(data.get("notes"), "Reason", NOTES_MAX_LENGTH)

    item = services.inventory_ledger.dispose_item(item.id, current_user.id, notes=notes)
    db.session.commit()
    return item


# --- Locations ---------------------------------------------------------------


def _read_location_fields(data: Dict[str, Any]) -> Dict[str, Optional[str]]:
    """Reads a location's name and description."""
    return {
        "name": _read_required_text(data.get("name"), "Name", LOCATION_NAME_MAX_LENGTH),
        "description": _read_text(
            data.get("description"), "Description", LOCATION_DESCRIPTION_MAX_LENGTH
        ),
    }


def _check_location_name_is_unique(
    workgroup: models.WorkGroup, name: str, exclude_id: Optional[int] = None
) -> None:
    """Location names must be unique within a workgroup, ignoring case."""
    query = db.session.query(models.Location).filter(
        models.Location.workgroup_id == workgroup.id,
        func.lower(models.Location.name) == name.lower(),
    )
    if exclude_id is not None:
        query = query.filter(models.Location.id != exclude_id)
    if query.first() is not None:
        raise InventoryValidationError("There's already a location with this name.")


def add_location(workgroup: models.WorkGroup, data: Dict[str, Any]) -> models.Location:
    """
    Adds a storage location.

    Args:
        workgroup: The workgroup.
        data: The JSON body, with `name` and optional `description`.

    Returns:
        The new location.
    """
    fields = _read_location_fields(data)
    _check_location_name_is_unique(workgroup, fields["name"])

    location = models.Location(workgroup_id=workgroup.id, is_archived=False, **fields)
    db.session.add(location)
    db.session.commit()
    return location


def edit_location(
    workgroup: models.WorkGroup, location_id: int, data: Dict[str, Any]
) -> models.Location:
    """
    Renames a location or changes its description.

    Args:
        workgroup: The workgroup.
        location_id: The Location id.
        data: The JSON body, with `name` and optional `description`.

    Returns:
        The updated location.
    """
    location = _get_location(workgroup, location_id)
    fields = _read_location_fields(data)
    _check_location_name_is_unique(workgroup, fields["name"], exclude_id=location.id)

    location.name = fields["name"]
    location.description = fields["description"]
    db.session.commit()
    return location


def set_location_archived(
    workgroup: models.WorkGroup, location_id: int, data: Dict[str, Any]
) -> models.Location:
    """
    Archives or restores a location. A location that still holds active
    containers can't be archived.

    Args:
        workgroup: The workgroup.
        location_id: The Location id.
        data: The JSON body, with `isArchived` as a boolean.

    Returns:
        The updated location.
    """
    location = _get_location(workgroup, location_id)
    is_archived = data.get("isArchived")
    if not isinstance(is_archived, bool):
        raise InventoryValidationError("isArchived must be true or false.")

    if is_archived:
        active_containers = location.inventory_items.filter(
            models.InventoryItem.is_archived.is_(False)
        ).count()
        if active_containers:
            raise InventoryValidationError(
                "Move or dispose of this location's containers before archiving it."
            )

    location.is_archived = is_archived
    db.session.commit()
    return location


# --- Compound lookup and novel compounds -------------------------------------


def _user_novel_compounds_query(workgroup: models.WorkGroup):
    """Novel compounds in the workbooks the user belongs to in this workgroup."""
    workbook_ids = [
        wb.id
        for wb in services.workbook.get_workbooks_from_user_group_combination(
            workgroup.name
        )
    ]
    return db.session.query(models.NovelCompound).filter(
        models.NovelCompound.workbook.in_(workbook_ids)
    )


def _suggestions(workgroup: models.WorkGroup, text: str) -> List[str]:
    """Up to MAX_SUGGESTIONS compound names similar to `text`, novel compounds first."""
    pattern = text.lower()
    # autoescape stops % and _ in the user's text acting as SQL wildcards
    novel_names = [
        nc.name
        for nc in _user_novel_compounds_query(workgroup)
        .filter(
            func.lower(models.NovelCompound.name).contains(pattern, autoescape=True)
        )
        .order_by(models.NovelCompound.name)
        .limit(MAX_SUGGESTIONS)
    ]
    compound_names = [
        compound.name
        for compound in db.session.query(models.Compound)
        .filter(func.lower(models.Compound.name).startswith(pattern, autoescape=True))
        .order_by(models.Compound.name)
        .limit(MAX_SUGGESTIONS)
    ]
    return list(dict.fromkeys(novel_names + compound_names))[:MAX_SUGGESTIONS]


def lookup_compound(
    workgroup: models.WorkGroup, query: Optional[str]
) -> Dict[str, Any]:
    """
    Looks a compound up by CAS number or name, in the Compound database and the
    novel compounds of the user's workbooks. Read-only: unlike the reaction
    table's lookup, it doesn't add to any recent compounds list.

    Args:
        workgroup: The workgroup.
        query: A CAS number or compound name.

    Returns:
        `{status: "found", compound}`, `{status: "cas_not_found", cas}` or
        `{status: "not_found", suggestions}`.
    """
    text = (query or "").strip()
    if not text:
        raise InventoryValidationError("Enter a name or CAS number.")
    if len(text) > MAX_LOOKUP_QUERY_LENGTH:
        raise InventoryValidationError("That search is too long.")

    novel_compounds = _user_novel_compounds_query(workgroup)
    if CAS_PATTERN.fullmatch(text):
        chemical = (
            services.compound.from_cas(text)
            or novel_compounds.filter(models.NovelCompound.cas == text).first()
        )
        if chemical is None:
            return {"status": "cas_not_found", "cas": text}
    else:
        chemical = (
            services.compound.from_name(text)
            or novel_compounds.filter(
                func.lower(models.NovelCompound.name) == text.lower()
            ).first()
        )
        if chemical is None:
            return {"status": "not_found", "suggestions": _suggestions(workgroup, text)}

    return {"status": "found", "compound": serialise_chemical(chemical, workgroup)}


def _form_number(value: Optional[float]) -> str:
    """Formats an optional number the way the novel compound form sends it."""
    return "" if value is None else repr(value)


def add_novel_compound(
    workgroup: models.WorkGroup, data: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Adds a novel compound to one of the user's workbooks, with the same
    validation as the reaction table (unique name, CAS and structure, valid
    SMILES and hazard codes).

    Args:
        workgroup: The workgroup.
        data: The JSON body: `workbookId`, `name`, and optional `cas`, `smiles`,
            `molWeight`, `density`, `concentration` and `hazardCodes`.

    Returns:
        The new compound, serialised as for the lookup.
    """
    workbook_id = _read_id(data.get("workbookId"), "Workbook")
    workbook = next(
        (
            wb
            for wb in services.workbook.get_workbooks_from_user_group_combination(
                workgroup.name
            )
            if wb.id == workbook_id
        ),
        None,
    )
    if workbook is None:
        raise InventoryNotFoundError()

    name = _read_required_text(data.get("name"), "Name", 200)
    cas = _read_text(data.get("cas"), "CAS number", 20)
    if cas and not CAS_PATTERN.fullmatch(cas):
        raise InventoryValidationError("Enter a CAS number like 67-64-1.")
    smiles = _read_text(data.get("smiles"), "SMILES", 1000)
    mol_weight = _read_optional_positive(
        data.get("molWeight"), "Molecular weight", 100000
    )
    density = _read_optional_positive(data.get("density"), "Density", 100)
    concentration = _read_optional_positive(
        data.get("concentration"), "Concentration", 1000
    )
    hazard_codes = _read_text(data.get("hazardCodes"), "Hazard codes", 500)

    new_compound = services.novel_compound.NewNovelCompound(
        workbook,
        form={
            "name": name,
            "cas": cas or "",
            "smiles": smiles or "",
            "molWeight": _form_number(mol_weight),
            "density": _form_number(density),
            "concentration": _form_number(concentration),
            "hPhrase": hazard_codes or "",
            # anything but "sketcher" means the SMILES is validated
            "source": "inventory",
        },
    )
    if new_compound.validate() != "validation_successful":
        raise InventoryValidationError(new_compound.feedback)

    novel_compound = services.novel_compound.add(
        new_compound.name,
        new_compound.cas,
        new_compound.mol_formula,
        new_compound.mol_weight,
        new_compound.density,
        new_compound.concentration,
        new_compound.hazard_codes,
        new_compound.smiles,
        new_compound.inchi,
        new_compound.inchi_key,
        workbook.id,
    )
    return serialise_chemical(novel_compound, workgroup)
