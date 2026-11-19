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

    A step records three kinds of value, all defined per step type in
    services.workup_steps.STEP_DEFINITIONS and validated before saving:
    - `materials`: what was used (solvents, drying agents, ...), as linked rows with
      amounts, so green chemistry metrics can query them.
    - `temperature`: in °C, a column of its own so it can be queried too.
    - `data`: the step's other settings (duration, number of washes, ...) as JSON.
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
    temperature = db.Column(db.Numeric(7, 2), nullable=True)
    data = db.Column(db.JSON, nullable=False, default=dict)
    notes = db.Column(db.Text, nullable=False, default="")
    time_of_creation = db.Column(db.DateTime, nullable=False, default=_now)
    time_of_update = db.Column(db.DateTime, nullable=False, default=_now, onupdate=_now)
    materials = db.relationship(
        "WorkupStepMaterial",
        order_by="[WorkupStepMaterial.slot, WorkupStepMaterial.sort_order]",
        back_populates="step",
        cascade="all, delete-orphan",
        passive_deletes=True,
    )

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
            "temperature": (
                float(self.temperature) if self.temperature is not None else None
            ),
            "data": self.data or {},
            "materials": [material.to_dict() for material in self.materials],
            "notes": self.notes or "",
        }
