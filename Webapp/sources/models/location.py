from sources.extensions import db
from sqlalchemy import UniqueConstraint


class Location(db.Model):
    """A physical storage location, scoped to a workgroup."""

    __tablename__ = "Location"
    __table_args__ = (
        UniqueConstraint("workgroup_id", "name", name="uq_location_workgroup_name"),
    )

    id = db.Column(db.Integer, primary_key=True)
    workgroup_id = db.Column(
        db.Integer, db.ForeignKey("WorkGroup.id"), nullable=False, index=True
    )
    name = db.Column(db.String(255), nullable=False)
    description = db.Column(db.Text)
    is_archived = db.Column(db.Boolean, default=False, nullable=False)

    # Later, for hierarchy (building > room > cabinet > shelf):
    # parent_id = db.Column(db.Integer, db.ForeignKey("location.id"), nullable=True)

    workgroup = db.relationship("WorkGroup", backref="locations")
    inventory_items = db.relationship(
        "InventoryItem", back_populates="location", lazy="dynamic"
    )
