from sources.extensions import db
from sources.models.data_export_request import ApprovalStatus
from sources.models.template_instance import InstanceType, TemplateInstance


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
