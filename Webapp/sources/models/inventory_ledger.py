import enum

from sources.extensions import db
from sqlalchemy import CheckConstraint, Index


class StockChangeReason(str, enum.Enum):
    RECEIVED = "received"
    USED_IN_REACTION = "used_in_reaction"
    ADJUSTMENT = "adjustment"
    DISPOSED = "disposed"


class InventoryLedgerEntry(db.Model):
    """Append-only record of every quantity change to an InventoryItem.

    Never update or delete rows. To correct a mistake, add a reversing entry.
    `change` is always in the item's unit.
    """

    __tablename__ = "InventoryLedgerEntry"
    __table_args__ = (
        CheckConstraint("change <> 0", name="ck_ledger_change_nonzero"),
        Index("ix_ledger_item_created", "inventory_item_id", "created_at"),
    )

    id = db.Column(db.Integer, primary_key=True)
    inventory_item_id = db.Column(
        db.Integer, db.ForeignKey("InventoryItem.id"), nullable=False
    )
    change = db.Column(db.Numeric(12, 4), nullable=False)  # signed

    reason = db.Column(
        db.Enum(
            StockChangeReason,
            name="stock_change_reason",
            values_callable=lambda e: [m.value for m in e],
        ),
        nullable=False,
    )

    reaction_id = db.Column(
        db.Integer, db.ForeignKey("Reaction.id"), nullable=True  # adjust table name
    )
    user_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    notes = db.Column(db.Text)
    created_at = db.Column(
        db.DateTime(timezone=True), server_default=db.func.now(), nullable=False
    )

    inventory_item = db.relationship("InventoryItem", back_populates="ledger_entries")
