from typing import Any, Dict, Optional

from sources.extensions import db
from sqlalchemy import CheckConstraint

from ..base import Model
from ..solvent import Solvent
from .enums import WorkupMaterialFate, WorkupMaterialRole


class WorkupStepMaterial(Model):
    """
    A material used in a workup step (a solvent, drying agent, stationary phase, ...),
    kept as its own row so amounts can be summed for green chemistry metrics.

    The material is linked to at most one of: the Compound database
    (`compound_id`), a workbook's novel compound (`novel_compound_id`), or a CHEM21
    solvent that has no Compound row (`solvent_name`). With none set it is unlinked
    and only `name` describes it (e.g. "brine").

    `quantity` is the amount used each time. Some slots are repeated by the step's
    `repeats` value (e.g. three washes), so totals must come from
    services.workup_steps.total_quantity rather than `quantity` alone.
    """

    __tablename__ = "WorkupStepMaterial"
    __table_args__ = (
        # an "exclusive arc": at most one link; none for unidentified materials
        CheckConstraint(
            "num_nonnulls(compound_id, novel_compound_id, solvent_name) <= 1",
            name="ck_workup_step_material_one_source",
        ),
    )

    id = db.Column(db.Integer, primary_key=True)
    step_id = db.Column(
        db.Integer,
        db.ForeignKey("WorkupStep.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    slot = db.Column(db.Text, nullable=False)
    sort_order = db.Column(db.Integer, nullable=False)
    role = db.Column(db.Enum(WorkupMaterialRole), nullable=False)
    fate = db.Column(
        db.Enum(WorkupMaterialFate),
        nullable=False,
        default=WorkupMaterialFate.WASTE,
        server_default=WorkupMaterialFate.WASTE.value,
    )
    compound_id = db.Column(
        db.Integer,
        db.ForeignKey("Compound.id", ondelete="SET NULL"),
        nullable=True,
        index=True,
    )
    novel_compound_id = db.Column(
        db.Integer,
        db.ForeignKey("NovelCompound.id", ondelete="SET NULL"),
        nullable=True,
        index=True,
    )
    solvent_name = db.Column(
        db.Text,
        db.ForeignKey("Solvent.name", ondelete="SET NULL", onupdate="CASCADE"),
        nullable=True,
        index=True,
    )
    name = db.Column(db.Text, nullable=False, default="")
    quantity = db.Column(db.Numeric(12, 4), nullable=True)
    unit = db.Column(db.String(16), nullable=False)
    # percentage of a mixture, e.g. 90 for the hexane in 90:10 hexane/EtOAc; only
    # for slots with proportions, where a step's values add up to 100
    proportion = db.Column(db.Numeric(5, 2), nullable=True)

    step = db.relationship("WorkupStep", back_populates="materials")
    compound = db.relationship("Compound")
    novel_compound = db.relationship("NovelCompound")
    solvent = db.relationship("Solvent")

    @property
    def chem21_flag(self) -> Optional[int]:
        """
        The CHEM21 sustainability flag (1 highly hazardous ... 4 recommended, 5 not
        in CHEM21) when the material is a known solvent, otherwise None.
        """
        if self.solvent is not None:
            return self.solvent.flag
        # Compound.solvent and NovelCompound.solvent hold the linked Solvent's name
        chemical = self.compound or self.novel_compound
        if chemical is None or not chemical.solvent:
            return None
        solvent = db.session.get(Solvent, chemical.solvent)
        return solvent.flag if solvent is not None else None

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialises the material for the workup card.

        Returns:
            Dictionary representation of the material.
        """
        return {
            "id": self.id,
            "slot": self.slot,
            "role": self.role.value,
            "fate": self.fate.value,
            "compoundId": self.compound_id,
            "novelCompoundId": self.novel_compound_id,
            "solventName": self.solvent_name,
            "name": self.name or "",
            "quantity": float(self.quantity) if self.quantity is not None else None,
            "unit": self.unit,
            "proportion": (
                float(self.proportion) if self.proportion is not None else None
            ),
            "flag": self.chem21_flag,
        }
