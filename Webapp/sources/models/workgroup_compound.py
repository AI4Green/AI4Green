from datetime import datetime, timezone

from sources.extensions import db
from sqlalchemy import UniqueConstraint


class WorkgroupCompound(db.Model):
    """A compound a workgroup holds (or tracks). One row per (workgroup, compound)."""

    __tablename__ = "WorkgroupCompound"
    __table_args__ = (
        UniqueConstraint("workgroup_id", "compound_id", name="uq_workgroup_compound"),
    )

    id = db.Column(db.Integer, primary_key=True)
    workgroup_id = db.Column(
        db.Integer, db.ForeignKey("WorkGroup.id"), nullable=False, index=True
    )
    compound_id = db.Column(
        db.Integer, db.ForeignKey("Compound.id"), nullable=False, index=True
    )

    # Workgroup-level info (true of the compound for the lab, not per bottle)
    preferred_supplier = db.Column(db.String(255))
    min_stock_threshold = db.Column(db.Numeric(12, 4))  # in threshold_unit
    threshold_unit = db.Column(db.String(16))
    notes = db.Column(db.Text)
    created_at = db.Column(
        db.DateTime(timezone=True), default=db.func.now(), nullable=False
    )

    compound = db.relationship("Compound", backref="workgroup_compounds")
    workgroup = db.relationship("WorkGroup", backref="workgroup_compounds")
    inventory_items = db.relationship(
        "InventoryItem", back_populates="workgroup_compound", lazy="dynamic"
    )
