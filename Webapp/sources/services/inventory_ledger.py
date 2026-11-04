from decimal import Decimal

from models.inventory_ledger import StockChangeReason
from sources import db, models


def adjust_stock(
    item_id,
    change,
    reason,
    user_id,
    reaction_id=None,
    notes=None,
):
    """The only place that should ever write InventoryItem.quantity.

    Locks the item row, checks the result isn't negative, writes the ledger
    entry and updates the cached quantity. The caller commits.
    """
    reason = StockChangeReason(reason)  # raises ValueError if invalid
    change = Decimal(str(change))
    if change == 0:
        raise ValueError("Change must be non-zero.")
    if reason is StockChangeReason.USED_IN_REACTION and reaction_id is None:
        raise ValueError("A reaction_id is required for used_in_reaction.")

    # Row lock so two concurrent changes to the same bottle queue up
    item = (
        db.session.query(models.InventoryItem)
        .filter_by(id=item_id)
        .with_for_update()
        .one()
    )

    new_quantity = item.quantity + change
    if new_quantity < 0:
        raise ValueError(
            f"Insufficient stock: have {item.quantity} {item.unit}, "
            f"change is {change}."
        )

    entry = models.InventoryLedgerEntry(
        inventory_item_id=item.id,
        change=change,
        reason=reason.value,
        reaction_id=reaction_id,
        user_id=user_id,
        notes=notes,
    )
    item.quantity = new_quantity
    db.session.add(entry)
    db.session.flush()
    return entry


def create_inventory_item(
    workgroup_compound_id, location_id, quantity, unit, user_id, **fields
):
    """Create a container with quantity 0, then record the initial amount as 'received'.

    `fields` takes the optional columns (supplier, lot_number, purity, barcode,
    date_received, date_opened, expiry_date).
    """
    item = models.InventoryItem(
        workgroup_compound_id=workgroup_compound_id,
        location_id=location_id,
        quantity=Decimal("0"),
        unit=unit,
        **fields,
    )
    if not item.belongs_to_single_workgroup():
        raise ValueError("Compound and location must belong to the same workgroup.")

    db.session.add(item)
    db.session.flush()  # get item.id
    adjust_stock(item.id, quantity, StockChangeReason.RECEIVED, user_id)
    return item


def dispose_item(item_id, user_id, notes=None):
    """Write off the remaining quantity and archive the container."""
    item = (
        db.session.query(models.InventoryItem)
        .filter_by(id=item_id)
        .with_for_update()
        .one()
    )
    if item.quantity > 0:
        adjust_stock(
            item.id, -item.quantity, StockChangeReason.DISPOSED, user_id, notes=notes
        )
    item.is_archived = True
    item.disposed_date = db.func.current_date()
    db.session.flush()
    return item
