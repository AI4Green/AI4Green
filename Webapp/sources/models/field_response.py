from sources.extensions import db

from .base import Model


class FieldResponse(Model):
    __tablename__ = "FieldResponse"

    id = db.Column(db.Integer, primary_key=True)

    field_id = db.Column(db.Integer, db.ForeignKey("Field.id"), nullable=False)
    field = db.relationship("Field", back_populates="field_response")

    field_response_values = db.relationship(
        "FieldResponseValue", back_populates="field_response"
    )

    comment = db.relationship("Comment", back_populates="field_response")

    template_instance_id = db.Column(
        db.Integer, db.ForeignKey("TemplateInstance.id"), nullable=False
    )
    template_instance = db.relationship(
        "TemplateInstance", back_populates="field_responses"
    )

    def to_dict(self):
        return {
            "id": self.id,
            "field_id": self.field_id,
            "field_response_values": self.field_response_values,  # does this get updates with multiple values?
            "template_instance_id": self.template_instance_id,
            "comment": self.comment,
        }
