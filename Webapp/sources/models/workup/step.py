from datetime import datetime
from typing import Any, Dict

import pytz
from sources.extensions import db

from ..base import Model
from .enums import WorkupCategory, WorkupStepType


def _now() -> datetime:
    """Returns the current London time without tzinfo, as the other models store it."""
    return datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)


class WorkupStep(Model):
    """
    One step of a reaction's workup or purification.

    The step's values live in `data`, a JSON object whose keys are fixed by the step
    type (see services.workup_steps.STEP_DEFINITIONS) and validated before saving.
    """

    __tablename__ = "WorkupStep"

    id = db.Column(db.Integer, primary_key=True)
    reaction_id = db.Column(
        db.Integer,
        db.ForeignKey("Reaction.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    reaction = db.relationship(
        "Reaction",
        backref=db.backref(
            "workup_steps",
            order_by="WorkupStep.sort_order",
            cascade="all, delete-orphan",
            passive_deletes=True,
        ),
    )
    category = db.Column(db.Enum(WorkupCategory), nullable=False)
    sort_order = db.Column(db.Integer, nullable=False)
    step_type = db.Column(db.Enum(WorkupStepType), nullable=False)
    data = db.Column(db.JSON, nullable=False, default=dict)
    notes = db.Column(db.Text, nullable=False, default="")
    time_of_creation = db.Column(db.DateTime, nullable=False, default=_now)
    time_of_update = db.Column(db.DateTime, nullable=False, default=_now, onupdate=_now)

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialises the step for the workup card.

        Returns:
            Dictionary representation of the step.
        """
        return {
            "id": self.id,
            "category": self.category.value,
            "sortOrder": self.sort_order,
            "stepType": self.step_type.value,
            "data": self.data or {},
            "notes": self.notes or "",
        }
