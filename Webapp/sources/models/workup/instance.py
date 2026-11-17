from typing import Any, Dict, List

from sources.extensions import db
from sources.models.template_instance import InstanceType, TemplateInstance


class WorkupInstance(TemplateInstance):
    """
    A workup template imported into a reaction.

    The instance only stores what it changes (see models.InstanceSectionOverride and
    Section.template_instance_id). Its steps are the template's steps, in the
    instance's order, plus any steps it added.
    """

    __tablename__ = "WorkupInstance"

    id = db.Column(
        db.Integer,
        db.ForeignKey("TemplateInstance.id", ondelete="CASCADE"),
        primary_key=True,
    )

    # steps this instance added for itself; the template is not changed
    extra_sections = db.relationship(
        "Section",
        back_populates="workup_instance",
        cascade="all, delete-orphan",
    )

    # this instance's position and hidden state for its steps; only written once
    # a user reorders or hides a step
    section_overrides = db.relationship(
        "InstanceSectionOverride",
        back_populates="template_instance",
        cascade="all, delete-orphan",
    )

    __mapper_args__ = {
        "polymorphic_identity": InstanceType.WORKUP,
    }

    def to_summary_dict(self) -> Dict[str, Any]:
        """
        Adds the template's category (workup or purification) to the summary.

        Returns:
            Summary representation of this instance.
        """
        data = super().to_summary_dict()
        data["category"] = self.template.category.value
        return data

    def ordered_sections(self) -> List[Dict[str, Any]]:
        """
        Merges the template's steps with this instance's changes.

        Each step uses the instance's override (position and hidden flag) if there
        is one, otherwise its own sort_order. Ties are broken by section id, so the
        order is stable.

        Returns:
            Serialised steps in display order. Each has "origin" ("template" or
            "instance") and "isHidden".
        """
        overrides = {x.section_id: x for x in self.section_overrides}
        sources = [(section, "template") for section in self.template.sections] + [
            (section, "instance") for section in self.extra_sections
        ]

        merged = []
        for section, origin in sources:
            override = overrides.get(section.id)
            data = section.to_dict(self.id)
            data["origin"] = origin
            data["isHidden"] = bool(override and override.is_hidden)
            if override is not None:
                data["sortOrder"] = override.sort_order
            merged.append(data)

        return sorted(merged, key=lambda x: (x["sortOrder"] or 0, x["id"]))

    def to_dict(self) -> Dict[str, Any]:
        """
        Full representation, with the merged step list in place of the template's sections.

        Returns:
            Full representation of this instance.
        """
        data = self.to_summary_dict()
        data.update(
            {
                "sections": self.ordered_sections(),
                "fieldResponses": [x.to_dict() for x in self.field_responses],
            }
        )
        return data
