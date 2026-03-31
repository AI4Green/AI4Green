from enum import Enum

from sources.extensions import db

from .base import Model


class TemplateType(Enum):
    """
    Enum for templating types, allows generic template db table
    """

    COSHH = "COSHH"


class ApprovalStatus(Enum):
    """
    Enum for Approval Status, todo can we reuse a previous one?
    """

    PENDING = "PENDING"
    APPROVED = "APPROVED"
    REJECTED = "REJECTED"
    CHANGES_SUGGESTED = "CHANGES_SUGGESTED"


class Template(Model):
    """
    For more specific templates, you can inherit from this class using alembic polymorphism:
    needs specific TemplateType Eum value. You will also need to uncomment out the polymorphic config in this model

    ps: I havent tested this
    eg/

    Class NewTemplate(Template):
        __mapper_args__ = {"polymorphic_identity": TemplateType.NewTemplate}
    """

    __tablename__ = "Template"
    # __bind_key__ = "main" do we need this to add to the default db?

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    description = db.Column(db.Text, nullable=False, default="")
    template_type = db.Column(db.Enum(TemplateType), nullable=False)
    time_of_creation = db.Column(db.Datetime, nullable=False)
    timE_of_update = db.Column(db.Datetime)

    institution_id = db.Column(
        db.Integer, db.ForeignKey("Institution.id", ondelete="CASCADE")
    )
    institution = db.relationship("Institution", backref="template")

    sections = db.relationship("Section", back_populates="template")
    instances = db.relationship("Instance", back_populates="template")

    approved = db.Column(db.Enum(ApprovalStatus), nullable=False)

    # # Polymorphic config, uncomment for additional child template
    # __mapper_args__ = {"polymorphic_on": template_type}
