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

    reaction_id = db.Column(db.Integer, db.ForeignKey("Reaction.id"), nullable=False)
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

    def to_dict(self):
        workbook = self.reaction.workbook
        return {
            "id": self.id,
            "uuid": self.uuid,
            # "template_type": self.template_type.value,
            "owner_id": self.owner_id,
            "template_id": self.template_id,
            "reaction_id": self.reaction_id,
            "approver_id": self.approver_id,
            "sections": [x.to_dict() for x in self.template.sections],
            "fieldResponses": [x.to_dict() for x in self.field_responses],
            "approvalStatus": self.approval_status.value[0],
            "reactionCode": self.reaction.reaction_id,
            "workbook": workbook.name,
            "workgroup": workbook.WorkGroup.name,
        }
