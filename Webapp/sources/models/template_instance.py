from enum import Enum

from sources.extensions import db

from .base import Model


class InstanceType(Enum):
    GENERIC = "GENERIC"
    COSHH = "COSHH"


class TemplateInstance(Model):
    __tablename__ = "TemplateInstance"

    id = db.Column(db.Integer, primary_key=True)
    uuid = db.Column(db.Text)  # identifier, needed? or just use the id?

    owner_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    owner = db.relationship(
        "User", backref="template_instances", foreign_keys=[owner_id]
    )  # backref used so we dont have to edit user table

    template_id = db.Column(db.Integer, db.ForeignKey("Template.id"), nullable=False)
    template = db.relationship("Template", back_populates="template_instances")

    reaction_id = db.Column(db.Integer, db.ForeignKey("Reaction.id"), nullable=True)
    reaction = db.relationship("Reaction", backref="template_instances")

    field_responses = db.relationship(
        "FieldResponse", back_populates="template_instance"
    )

    instance_type = db.Column(
        db.Enum(InstanceType),
        nullable=False,
    )
    # Polymorphic config, uncomment for additional child template
    __mapper_args__ = {
        "polymorphic_on": instance_type,
        "polymorphic_identity": InstanceType.GENERIC,
    }

    def to_summary_dict(self):
        """
        Lightweight dictionary representation, intended for list/search views.

        Contains only the fields common to every instance type (identity, ownership,
        template and reaction/workbook/workgroup context). Does NOT include sections
        or field responses, to keep list payloads small and avoid extra queries.

        Inheritance:
            Subclasses should override this method, call ``super().to_summary_dict()``
            and extend the returned dict with their own summary fields. Because
            ``to_dict()`` is built on top of this method, any fields added here by a
            subclass automatically appear in the full representation too.

        Returns:
            dict: summary representation of this instance.
        """
        reaction_code = None
        workbook_name = None
        workgroup_name = None

        if self.reaction_id:
            reaction_code = self.reaction.reaction_id
            workbook = self.reaction.workbook
            workbook_name = workbook.name
            workgroup_name = workbook.WorkGroup.name

        return {
            "id": self.id,
            "uuid": self.uuid,
            "instance_type": self.instance_type.value,
            "owner_id": self.owner_id,
            "template_id": self.template_id,
            "template_name": self.template.name,
            "reaction_id": self.reaction_id,
            "reaction_code": reaction_code,
            "workbook": workbook_name,
            "workgroup": workgroup_name,
        }

    def to_dict(self):
        """
        Full dictionary representation, intended for detail views (a single form).

        Built as ``self.to_summary_dict()`` plus the heavy fields (template sections
        and field responses), so the result is always a superset of the summary.

        Inheritance:
            Subclasses should NOT normally override this method. ``self`` is the
            concrete subclass instance, so ``self.to_summary_dict()`` resolves to the
            subclass's override and its extra summary fields are included for free.
            Override only if the subclass needs extra heavy fields, and in that case
            call ``super().to_dict()`` and extend the result.

            Dependency direction must stay one-way (``to_dict`` -> ``to_summary_dict``).

        Returns:
            dict: full representation of this instance.
        """
        data = self.to_summary_dict()
        data.update(
            {
                "sections": [x.to_dict() for x in self.template.sections],
                "fieldResponses": [x.to_dict() for x in self.field_responses],
            }
        )
        return data
