from typing import Any, Dict

from sources.extensions import db

from ..base import Model
from .enums import AdditionDevice, AdditionForm, AdditionSpeed, MethodComponentRole
from .method import enum_value, number_value


class ReactionAddition(Model):
    """
    One addition to the reaction vessel: an ORD ReactionInput.

    Additions are numbered by `sort_order`. An addition with `with_previous` set
    is made at the same time as the one before it, so both share an ORD
    addition_order. `hold_time_min` is how long the mixture is left after this
    addition before the next one ("...and stirred for 30 min"); summing these
    gives each addition's ORD addition_time.
    """

    __tablename__ = "ReactionAddition"

    id = db.Column(db.Integer, primary_key=True)
    method_id = db.Column(
        db.Integer,
        db.ForeignKey("ReactionMethod.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    method = db.relationship("ReactionMethod", back_populates="additions")
    sort_order = db.Column(db.Integer, nullable=False)
    with_previous = db.Column(db.Boolean, nullable=False, default=False)
    form = db.Column(db.Enum(AdditionForm), nullable=True)
    speed = db.Column(db.Enum(AdditionSpeed), nullable=True)
    device = db.Column(db.Enum(AdditionDevice), nullable=True)
    duration_min = db.Column(db.Numeric(10, 2), nullable=True)
    temperature_c = db.Column(db.Numeric(7, 2), nullable=True)
    hold_time_min = db.Column(db.Numeric(10, 2), nullable=True)
    # only for solutions: molar concentration and the solvent it is dissolved in
    concentration_m = db.Column(db.Numeric(10, 4), nullable=True)
    solution_solvent_key = db.Column(db.Text, nullable=True)
    solution_solvent_name = db.Column(db.Text, nullable=True)
    notes = db.Column(db.Text, nullable=False, default="")

    components = db.relationship(
        "ReactionAdditionComponent",
        back_populates="addition",
        order_by="ReactionAdditionComponent.sort_order",
        cascade="all, delete-orphan",
        passive_deletes=True,
    )

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialises the addition with references to its components. Amounts are
        added by the service from the reaction table.

        Returns:
            Dictionary representation of the addition.
        """
        return {
            "id": self.id,
            "sortOrder": self.sort_order,
            "withPrevious": bool(self.with_previous),
            "form": enum_value(self.form),
            "speed": enum_value(self.speed),
            "device": enum_value(self.device),
            "durationMin": number_value(self.duration_min),
            "temperatureC": number_value(self.temperature_c),
            "holdTimeMin": number_value(self.hold_time_min),
            "concentrationM": number_value(self.concentration_m),
            "solutionSolventKey": self.solution_solvent_key,
            "solutionSolventName": self.solution_solvent_name,
            "notes": self.notes or "",
            "components": [component.to_dict() for component in self.components],
        }


class ReactionAdditionComponent(Model):
    """
    A reaction table component that is part of an addition. An addition has more
    than one when they are mixed before being added.

    The component is found in the reaction table by role and key (SMILES for
    reactants and reagents, the solvent id for solvents). `name` is kept so the
    method still reads sensibly if the component is later removed from the table.
    """

    __tablename__ = "ReactionAdditionComponent"

    id = db.Column(db.Integer, primary_key=True)
    addition_id = db.Column(
        db.Integer,
        db.ForeignKey("ReactionAddition.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    addition = db.relationship("ReactionAddition", back_populates="components")
    sort_order = db.Column(db.Integer, nullable=False)
    role = db.Column(db.Enum(MethodComponentRole), nullable=False)
    component_key = db.Column(db.Text, nullable=False)
    name = db.Column(db.Text, nullable=False, default="")

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialises the component reference.

        Returns:
            Dictionary representation of the component reference.
        """
        return {
            "id": self.id,
            "role": self.role.value,
            "key": self.component_key,
            "name": self.name,
        }
