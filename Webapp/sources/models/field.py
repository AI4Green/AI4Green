from sources.extensions import db

from .base import Model


class Field(Model):
    __tablename__ = "Field"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)
    sort_order = db.Column(db.Integer)
    mandatory = db.Column(db.Boolean, server_default=True)

    section_id = db.Column(db.Integer, db.ForeignKey("Section.id"))
    section = db.relationship("Section", back_populates="fields")

    field_response = db.relationship("FieldResponse", back_populates="field")

    select_field_options = db.relationship("SelectFieldOptions", back_populates="field")

    input_type_id = db.Column(db.Integer, db.ForeignKey("InputType.id"))
    input_type = db.relationship("InputType", back_populates="field")
