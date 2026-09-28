from sources.extensions import db
from sources.models.template import Template, TemplateType


class COSHHTemplate(Template):
    """
    Child template inherited from Template

    """

    __tablename__ = "COSHHTemplate"

    id = db.Column(
        db.Integer,
        db.ForeignKey("Template.id", ondelete="CASCADE"),
        primary_key=True,
    )

    sections = db.relationship(
        "Section",
        back_populates="template",
        cascade="all, delete-orphan",
    )

    # map back to parent template
    __mapper_args__ = {
        "polymorphic_identity": TemplateType.COSHH,
    }
