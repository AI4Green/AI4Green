from sources.extensions import db


class InventoryItem(db.Model):
    """A physical container of a compound, in a location."""

    __tablename__ = "InventoryItem"

    id = db.Column(db.Integer, primary_key=True)
    workgroup_compound_id = db.Column(
        db.Integer, db.ForeignKey("WorkgroupCompound.id"), nullable=False, index=True
    )
    workgroup_compound = db.relationship(
        "WorkgroupCompound", back_populates="inventory_items"
    )

    location_id = db.Column(
        db.Integer, db.ForeignKey("Location.id"), nullable=False, index=True
    )
    location = db.relationship("Location", back_populates="inventory_items")

    # Amount
    quantity = db.Column(db.Numeric(12, 4), nullable=False)
    unit = db.Column(db.String(16), nullable=False)  # e.g. g, mg, mL, L

    # Container details
    supplier = db.Column(db.String(255))
    lot_number = db.Column(db.String(100))
    purity = db.Column(db.String(50))
    barcode = db.Column(db.String(100), index=True)

    # Dates
    date_received = db.Column(db.Date)
    date_opened = db.Column(db.Date)
    expiry_date = db.Column(db.Date)

    is_archived = db.Column(db.Boolean, default=False, nullable=False)
    disposed_date = db.Column(db.Date)

    created_at = db.Column(
        db.DateTime(timezone=True), default=db.func.now(), nullable=False
    )
    updated_at = db.Column(
        db.DateTime(timezone=True),
        default=db.func.now(),
        onupdate=db.func.now(),
        nullable=False,
    )

    ledger_entries = db.relationship(
        "InventoryLedgerEntry",
        back_populates="inventory_item",
        order_by="InventoryLedgerEntry.created_at",
        lazy="dynamic",
    )

    def belongs_to_single_workgroup(self) -> bool:
        """Service-layer check: compound and location must share a workgroup."""
        return self.workgroup_compound.workgroup_id == self.location.workgroup_id
