from enum import Enum

from sources.extensions import db

from .base import Model


class TemplateType(Enum):
    """
    Enum for templating types, allows generic template db table
    """

    COSHH = "COSHH"
    REACTIONS = "REACTIONS"


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

    creator_id = db.Column(
        db.Integer, db.ForeignKey("User.id", ondelete="SET NULL"), nullable=True
    )
    creator = db.relationship("User", backref="templates")

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

    def to_dict(self):
        """
        Simple method to serialise template. Does not include sections, these should be retrieved by api call
        """
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "time_of_creation": self.time_of_creation,
            "time_of_update": self.time_of_update,
            "creator_id": self.creator_id,
            "institution_id": self.institution_id,
            "stage": self.status.value.capitalize(),
            "inUseCount": len(
                self.template_instances
            ),  # camel case for returning to front end
            "permissions": [
                "CanPublish"
            ],  # included here for dev todo: implement permission at admin level
        }
