from enum import Enum
from typing import Any, Dict

from sources.extensions import db
from sources.models.template import Template, TemplateType


class WorkupCategory(Enum):
    """
    Which reaction card a workup template is shown in.

    ORD has a single ReactionWorkup type covering both, so this only splits them
    for users; the template type stays TemplateType.WORKUP.
    """

    WORKUP = "WORKUP"
    PURIFICATION = "PURIFICATION"


class WorkupTemplate(Template):
    """
    A workup or purification template, shared within one workbook.
    """

    __tablename__ = "WorkupTemplate"

    id = db.Column(
        db.Integer,
        db.ForeignKey("Template.id", ondelete="CASCADE"),
        primary_key=True,
    )

    workbook_id = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False
    )
    workbook = db.relationship("WorkBook", backref="workup_templates")

    category = db.Column(
        db.Enum(WorkupCategory), nullable=False, default=WorkupCategory.WORKUP
    )

    # map back to parent template
    __mapper_args__ = {
        "polymorphic_identity": TemplateType.WORKUP,
    }

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialises the template with its workbook and category.

        Returns:
            Dictionary representation of the template.
        """
        data = super().to_dict()
        data.update(
            {
                "workbookId": self.workbook_id,
                "category": self.category.value,
            }
        )
        return data
