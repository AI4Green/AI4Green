from sources.extensions import db
from sources.models.template import Template, TemplateType


class WorkupTemplate(Template):
    """
    Child template inherited from Template
    """

    __tablename__ = "WorkupTemplate"

    id = db.Column(
        db.Integer,
        db.ForeignKey("Template.id", ondelete="CASCADE"),
        primary_key=True,
    )

    steps = db.relationship(
        "WorkupSection",
        back_populates="template",
        cascade="all, delete-orphan",
    )

    # map back to parent template
    __mapper_args__ = {
        "polymorphic_identity": TemplateType.WORKUP,
    }
