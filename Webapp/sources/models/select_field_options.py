from sources.extensions import db

from .base import Model


class SelectFieldOptions(Model):
    __tablename__ = "SelectFieldOptions"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)

    field_id = db.Column(db.Integer, db.ForeignKey("Field.id"), nullable=False)
    field = db.relationship("Field", back_populates="select_field_options")
