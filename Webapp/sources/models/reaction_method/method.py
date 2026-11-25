from datetime import datetime
from decimal import Decimal
from enum import Enum
from typing import Any, Dict, Optional

import pytz
from sources.extensions import db

from ..base import Model
from .enums import (
    Atmosphere,
    EnvironmentType,
    PressureControl,
    StirringMethod,
    StirringRate,
    TemperatureControl,
    VesselMaterial,
    VesselType,
)


def _now() -> datetime:
    """Returns the current London time without tzinfo, as the other models store it."""
    return datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)


def enum_value(value: Optional[Enum]) -> Optional[str]:
    """
    Serialises an optional enum column.

    Args:
        value: The enum member, or None when unset.

    Returns:
        The member's value, or None.
    """
    return value.value if value is not None else None


def number_value(value: Optional[Decimal]) -> Optional[float]:
    """
    Serialises an optional Numeric column, which SQLAlchemy returns as a Decimal.

    Args:
        value: The column value, or None when unset.

    Returns:
        The value as a float, or None.
    """
    return float(value) if value is not None else None


class ReactionMethod(Model):
    """
    How a reaction is run: the ORD ReactionSetup and ReactionConditions, plus the
    ordered additions (ORD ReactionInputs) in `additions`.

    There is at most one method per reaction. Amounts are not stored here: each
    addition points at reaction table components, whose masses, moles and volumes
    are read from the reaction table whenever the method is shown.

    Units are fixed so the values can be compared and queried: mL, °C, bar, rpm
    and minutes.
    """

    __tablename__ = "ReactionMethod"

    id = db.Column(db.Integer, primary_key=True)
    reaction_id = db.Column(
        db.Integer,
        db.ForeignKey("Reaction.id", ondelete="CASCADE"),
        nullable=False,
        unique=True,
    )
    reaction = db.relationship(
        "Reaction",
        backref=db.backref(
            "method",
            uselist=False,
            cascade="all, delete-orphan",
            passive_deletes=True,
        ),
    )

    # ORD ReactionSetup
    vessel_type = db.Column(db.Enum(VesselType), nullable=True)
    vessel_material = db.Column(db.Enum(VesselMaterial), nullable=True)
    vessel_volume_ml = db.Column(db.Numeric(10, 2), nullable=True)
    # lists of VesselPreparation / VesselAttachment values, checked by the service
    vessel_preparations = db.Column(db.JSON, nullable=False, default=list)
    vessel_attachments = db.Column(db.JSON, nullable=False, default=list)
    environment = db.Column(db.Enum(EnvironmentType), nullable=True)

    # ORD ReactionConditions
    temperature_control = db.Column(db.Enum(TemperatureControl), nullable=True)
    temperature_c = db.Column(db.Numeric(7, 2), nullable=True)
    pressure_control = db.Column(db.Enum(PressureControl), nullable=True)
    atmosphere = db.Column(db.Enum(Atmosphere), nullable=True)
    pressure_bar = db.Column(db.Numeric(8, 3), nullable=True)
    stirring_method = db.Column(db.Enum(StirringMethod), nullable=True)
    stirring_rate = db.Column(db.Enum(StirringRate), nullable=True)
    stirring_rpm = db.Column(db.Integer, nullable=True)
    reflux = db.Column(db.Boolean, nullable=False, default=False)
    reaction_time_min = db.Column(db.Numeric(10, 2), nullable=True)
    details = db.Column(db.Text, nullable=False, default="")

    time_of_creation = db.Column(db.DateTime, nullable=False, default=_now)
    time_of_update = db.Column(db.DateTime, nullable=False, default=_now, onupdate=_now)

    additions = db.relationship(
        "ReactionAddition",
        back_populates="method",
        order_by="ReactionAddition.sort_order",
        cascade="all, delete-orphan",
        passive_deletes=True,
    )

    def setup_dict(self) -> Dict[str, Any]:
        """
        Serialises the ORD ReactionSetup fields.

        Returns:
            The setup, with camelCase keys.
        """
        return {
            "vesselType": enum_value(self.vessel_type),
            "vesselMaterial": enum_value(self.vessel_material),
            "vesselVolumeMl": number_value(self.vessel_volume_ml),
            "vesselPreparations": list(self.vessel_preparations or []),
            "vesselAttachments": list(self.vessel_attachments or []),
            "environment": enum_value(self.environment),
        }

    def conditions_dict(self) -> Dict[str, Any]:
        """
        Serialises the ORD ReactionConditions fields.

        Returns:
            The conditions, with camelCase keys.
        """
        return {
            "temperatureControl": enum_value(self.temperature_control),
            "temperatureC": number_value(self.temperature_c),
            "pressureControl": enum_value(self.pressure_control),
            "atmosphere": enum_value(self.atmosphere),
            "pressureBar": number_value(self.pressure_bar),
            "stirringMethod": enum_value(self.stirring_method),
            "stirringRate": enum_value(self.stirring_rate),
            "stirringRpm": self.stirring_rpm,
            "reflux": bool(self.reflux),
            "reactionTimeMin": number_value(self.reaction_time_min),
            "details": self.details or "",
        }

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialises the setup and conditions. Additions are serialised separately
        because they need the reaction table to show their amounts.

        Returns:
            Dictionary representation of the method.
        """
        return {
            "id": self.id,
            "setup": self.setup_dict(),
            "conditions": self.conditions_dict(),
        }
