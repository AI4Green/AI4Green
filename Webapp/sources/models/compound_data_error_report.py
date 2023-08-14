from sources.extensions import db
from .base import Model


class CompoundDataErrorReport(Model):
    __tablename__ = "CompoundDataErrorReport"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    compound_name = db.Column(db.Text, nullable=False)
    compound = db.Column(
        db.ForeignKey("Compound.id", ondelete="CASCADE"), nullable=False, index=True
    )
    error_type = db.Column(db.Text, nullable=False)
    additional_info = db.Column(db.Text, nullable=False)
    time = db.Column(db.DateTime, nullable=False)
