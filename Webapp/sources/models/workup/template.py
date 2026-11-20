from typing import Any, Dict

from sources.extensions import db

from ..base import Model
from .enums import WorkupCategory, WorkupStepType
from .step import _now


class WorkupTemplate(Model):
    """
    A saved sequence of workup or purification steps, shared within one workbook.

    Only the step types and their order are saved; a reaction that starts from the
    template gets its own empty copies of the steps, with no link back.
    """

    __tablename__ = "WorkupTemplate"

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    description = db.Column(db.Text, nullable=False, default="")
    workbook_id = db.Column(
        db.Integer, db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False
    )
    workbook = db.relationship(
        "WorkBook",
        backref=db.backref(
            "workup_templates", cascade="all, delete-orphan", passive_deletes=True
        ),
    )
    category = db.Column(db.Enum(WorkupCategory), nullable=False)
    creator_id = db.Column(
        db.Integer, db.ForeignKey("User.id", ondelete="SET NULL"), nullable=True
    )
    creator = db.relationship("User")
    time_of_creation = db.Column(db.DateTime, nullable=False, default=_now)

    steps = db.relationship(
        "WorkupTemplateStep",
        back_populates="template",
        order_by="WorkupTemplateStep.sort_order",
        cascade="all, delete-orphan",
        passive_deletes=True,
    )

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialises the template and its step types.

        Returns:
            Dictionary representation of the template.
        """
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "category": self.category.value,
            "creatorId": self.creator_id,
            "creator": self.creator.fullname if self.creator else None,
            "steps": [step.step_type.value for step in self.steps],
        }


class WorkupTemplateStep(Model):
    """One step type, in order, within a workup template."""

    __tablename__ = "WorkupTemplateStep"

    id = db.Column(db.Integer, primary_key=True)
    template_id = db.Column(
        db.Integer,
        db.ForeignKey("WorkupTemplate.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    template = db.relationship("WorkupTemplate", back_populates="steps")
    sort_order = db.Column(db.Integer, nullable=False)
    step_type = db.Column(db.Enum(WorkupStepType), nullable=False)
