from sources.extensions import db

from .base import Model


class FieldResponseValue(Model):
    __tablename__ = "FieldResponseValue"

    id = db.Column(db.Integer, primary_key=True)
    value = db.Column(db.String, nullable=False)
    time_of_response = db.Column(db.DateTime, nullable=False)

    field_response_id = db.Column(db.Integer, db.ForeignKey("FieldResponse.id"))
    field_response = db.relationship(
        "FieldResponse", back_populates="field_response_values"
    )
