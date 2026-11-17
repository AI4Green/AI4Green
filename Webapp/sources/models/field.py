from typing import Any, Dict, Optional

from sources.extensions import db

from .base import Model


class Field(Model):
    """
    A field in a section.

    template_instance_id is only set on a field that one workup instance added to a
    template section. Template fields leave it empty and are shared by every instance.
    """

    __tablename__ = "Field"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)
    sort_order = db.Column(db.Integer)
    mandatory = db.Column(db.Boolean, default="true")

    section_id = db.Column(db.Integer, db.ForeignKey("Section.id"))
    section = db.relationship("Section", back_populates="fields")

    template_instance_id = db.Column(
        db.Integer, db.ForeignKey("WorkupInstance.id", ondelete="CASCADE")
    )

    field_response = db.relationship("FieldResponse", back_populates="field")

    select_field_options = db.relationship("SelectFieldOptions", back_populates="field")

    input_type_id = db.Column(db.Integer, db.ForeignKey("InputType.id"))
    input_type = db.relationship("InputType", back_populates="field")

    def to_dict(self, template_instance_id: Optional[int] = None) -> Dict[str, Any]:
        """
        Serialises the field.

        Args:
            template_instance_id: If given, include only this instance's responses.

        Returns:
            Dictionary representation of the field.
        """
        return {
            "id": self.id,
            "name": self.name,
            "sortOrder": self.sort_order,
            "mandatory": self.mandatory,
            "section_id": self.section_id,
            "templateInstanceId": self.template_instance_id,
            "inputType": self.input_type.to_dict(),
            "fieldResponse": [
                x.to_dict()
                for x in self.field_response
                if x.template_instance_id == template_instance_id
            ],
        }
