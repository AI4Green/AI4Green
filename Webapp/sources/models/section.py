from typing import Any, Dict, Optional

from sources.extensions import db

from .base import Model


class Section(Model):
    """A section of a COSHH template, holding an ordered list of fields."""

    __tablename__ = "Section"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)

    sort_order = db.Column(db.Integer)

    template_id = db.Column(db.Integer, db.ForeignKey("COSHHTemplate.id"))
    template = db.relationship("COSHHTemplate", back_populates="sections")

    fields = db.relationship(
        "Field", back_populates="section", cascade="all, delete-orphan"
    )

    def to_dict(self, template_instance_id: Optional[int] = None) -> Dict[str, Any]:
        """
        Serialises the section and its fields.

        Args:
            template_instance_id: If given, include that instance's field responses.

        Returns:
            Dictionary representation of the section.
        """
        return {
            "id": self.id,
            "name": self.name,
            "sortOrder": self.sort_order,
            "templateId": self.template_id,
            "fields": [
                x.to_dict(template_instance_id)
                for x in sorted(self.fields, key=lambda field: field.sort_order or 0)
            ],
        }
