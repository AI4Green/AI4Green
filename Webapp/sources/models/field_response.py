from sources.extensions import db

from .base import Model


class FieldResponse(Model):
    __tablename__ = "FieldResponse"

    id = db.Column(db.Integer, primary_key=True)

    field_id = db.Column(db.Integer, db.ForeignKey("Field.id"), nullable=False)
    field = db.relationship("Field", backref="field_response")

    field_response_values = db.relationship(
        "FieldResponseValue", back_populates="field_response"
    )

    comment = db.relationship("Comment", back_populates="field_response")

    instance_id = db.Column(db.Integer, db.ForeignKey("Instance.id"), nullable=False)
    instance = db.relationship("Instance", back_populates="field_responses")
