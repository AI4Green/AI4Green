from sources.extensions import db
from sqlalchemy import CheckConstraint, UniqueConstraint


class WorkgroupCompound(db.Model):
    """A compound a workgroup holds (or tracks). One row per (workgroup, compound).

    The compound is either from the Compound database (`compound_id`) or a
    novel compound from one of the workgroup's workbooks (`novel_compound_id`).
    Exactly one of the two is set.
    """

    __tablename__ = "WorkgroupCompound"
    __table_args__ = (
        UniqueConstraint("workgroup_id", "compound_id", name="uq_workgroup_compound"),
        UniqueConstraint(
            "workgroup_id", "novel_compound_id", name="uq_workgroup_novel_compound"
        ),
        CheckConstraint(
            "(compound_id IS NULL) <> (novel_compound_id IS NULL)",
            name="ck_workgroup_compound_one_source",
        ),
    )

    id = db.Column(db.Integer, primary_key=True)
    workgroup_id = db.Column(
        db.Integer, db.ForeignKey("WorkGroup.id"), nullable=False, index=True
    )
    compound_id = db.Column(
        db.Integer, db.ForeignKey("Compound.id"), nullable=True, index=True
    )
    novel_compound_id = db.Column(
        db.Integer, db.ForeignKey("NovelCompound.id"), nullable=True, index=True
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
    novel_compound = db.relationship("NovelCompound", backref="workgroup_compounds")
    workgroup = db.relationship("WorkGroup", backref="workgroup_compounds")
    inventory_items = db.relationship(
        "InventoryItem", back_populates="workgroup_compound", lazy="dynamic"
    )

    @property
    def chemical(self):
        """The Compound or NovelCompound this row tracks."""
        return self.compound or self.novel_compound
