from models.data_export_request import ApprovalStatus
from models.template_instance import InstanceType, TemplateInstance
from sources.extensions import db


class COSHHInstance(TemplateInstance):
    __tablename__ = "COSHHInstance"

    # approval
    approval_status = db.Column(
        db.Enum(ApprovalStatus), nullable=False, default=ApprovalStatus.DRAFT
    )
    approver_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
    approver = db.relationship(
        "User", backref="template_approvals", foreign_keys=[approver_id]
    )  # backref used so we dont have to edit user table

    __mapper_args__ = {
        "polymorphic_identity": InstanceType.COSHH,
    }
