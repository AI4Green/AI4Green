from typing import Any, Dict, List, Optional

from sources.extensions import db

from .base import Model


class Section(Model):
    """
    A section of a template, or a section added to a single template instance.

    A section belongs to exactly one owner:
    - template_id: a section defined on a template, shared by every instance of it.
    - template_instance_id: a section a workup instance added for itself (an extra
      step for one reaction). The template is never changed by this. Only workup
      instances can own sections, so this points at WorkupInstance.

    Polymorphic on section_type, so child section types (e.g. WorkupSection) can
    add their own columns.
    """

    __tablename__ = "Section"
    __table_args__ = (
        db.CheckConstraint(
            "(template_id IS NULL) <> (template_instance_id IS NULL)",
            name="section_has_one_owner",
        ),
    )

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String, nullable=False)

    sort_order = db.Column(db.Integer)

    template_id = db.Column(db.Integer, db.ForeignKey("Template.id"))
    template = db.relationship("Template", back_populates="sections")

    # named template_instance_id (points to workup instance) to match
    # FieldResponse.template_instance_id, which services compare against
    template_instance_id = db.Column(
        db.Integer, db.ForeignKey("WorkupInstance.id", ondelete="CASCADE")
    )
    workup_instance = db.relationship("WorkupInstance", back_populates="extra_sections")

    fields = db.relationship(
        "Field", back_populates="section", cascade="all, delete-orphan"
    )

    section_type = db.Column(db.String(50), nullable=False)

    __mapper_args__ = {
        "polymorphic_on": section_type,
        "polymorphic_identity": "section",
    }

    def fields_for_instance(self, template_instance_id: Optional[int]) -> List[Any]:
        """
        Gets the fields to show for one instance, in sort order.

        Template fields are shared by every instance. Fields an instance added
        to this section (Field.template_instance_id) are only shown to that instance.

        Args:
            template_instance_id: The instance being shown, or None for the template itself.

        Returns:
            The visible fields, sorted by sort_order.
        """
        visible = [
            field
            for field in self.fields
            if field.template_instance_id is None  # fields belong to the template
            or field.template_instance_id
            == template_instance_id  # fields belong to the instance
        ]
        return sorted(visible, key=lambda field: field.sort_order or 0)

    def to_dict(self, template_instance_id: Optional[int] = None) -> Dict[str, Any]:
        """
        Serialises the section and its fields.

        Args:
            template_instance_id: If given, include that instance's field responses
                and any fields the instance added.

        Returns:
            Dictionary representation of the section.
        """
        return {
            "id": self.id,
            "name": self.name,
            "sortOrder": self.sort_order,
            "templateId": self.template_id,
            "templateInstanceId": self.template_instance_id,
            "fields": [
                x.to_dict(template_instance_id)
                for x in self.fields_for_instance(template_instance_id)
            ],
        }
