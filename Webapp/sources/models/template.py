from enum import Enum

from sources.extensions import db

from .base import Model


class TemplateType(Enum):
    """
    Enum for templating types, allows generic template db table
    """

    COSHH = "COSHH"


class TemplateStatus(Enum):
    """
    Enum for publish status
    """

    PUBLISHED = "PUBLISHED"
    DRAFT = "DRAFT"


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

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    description = db.Column(db.Text, nullable=False, default="")
    template_type = db.Column(db.Enum(TemplateType), nullable=False)
    time_of_creation = db.Column(db.DateTime, nullable=False)
    time_of_update = db.Column(db.DateTime)

    institution_id = db.Column(
        db.Integer, db.ForeignKey("Institution.id", ondelete="CASCADE")
    )
    institution = db.relationship("Institution", backref="template")

    sections = db.relationship("Section", back_populates="template")
    template_instances = db.relationship("TemplateInstance", back_populates="template")

    status = db.Column(
        db.Enum(TemplateStatus), nullable=False, default=TemplateStatus.DRAFT
    )

    # # Polymorphic config, uncomment for additional child template
    # __mapper_args__ = {"polymorphic_on": template_type}
