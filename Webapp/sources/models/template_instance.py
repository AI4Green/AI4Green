from enum import Enum

from sources.extensions import db

from .base import Model


class InstanceType(Enum):
    COSHH = "COSHH"
    REACTION = "REACTION"  # included here for future use


class ApprovalStatus(Enum):
    """
    Enum for Approval Status, todo can we reuse a previous one?
    """

    DRAFT = "DRAFT"
    PENDING = "PENDING"
    APPROVED = "APPROVED"
    REJECTED = "REJECTED"
    CHANGES_SUGGESTED = "CHANGES_SUGGESTED"


class TemplateInstance(Model):
    __tablename__ = "TemplateInstance"

    id = db.Column(db.Integer, primary_key=True)
    uuid = db.Column(db.Text)  # identifier, needed? or just use the id?
    template_type = db.Columns(db.Enum(InstanceType))

    owner_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    owner = db.relationship(
        "User", backref="template_instances"
    )  # backref used so we dont have to edit user table

    template_id = db.Column(db.Integer, db.ForeignKey("Template.id"), nullable=False)
    template = db.relationship("Template", back_populates="template_instances")

    reaction_id = db.Column(db.Integer, db.ForeignKey("Reaction.id"), nullable=False)
    reaction = db.relationship("Reaction", back_populates="template_instances")

    # approval
    approval_status = db.Column(
        db.Enum(ApprovalStatus), nullable=False, deafult=ApprovalStatus.DRAFT
    )
    approver_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    approver = db.relationship(
        "User", backref="template_approvals"
    )  # backref used so we dont have to edit user table
