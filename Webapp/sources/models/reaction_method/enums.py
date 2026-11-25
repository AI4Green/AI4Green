"""
Enums for reaction methods.

Most follow the Open Reaction Database (ORD) schema, keeping its value names so a
method can later be exported as ORD ReactionSetup, ReactionConditions and
ReactionInput messages. ORD's UNSPECIFIED values are not repeated here: an unset
column (None) means unspecified.
"""

from enum import Enum


class VesselType(Enum):
    """ORD Vessel.VesselType: what the reaction is run in."""

    CUSTOM = "CUSTOM"
    ROUND_BOTTOM_FLASK = "ROUND_BOTTOM_FLASK"
    VIAL = "VIAL"
    WELL_PLATE = "WELL_PLATE"
    MICROWAVE_VIAL = "MICROWAVE_VIAL"
    TUBE = "TUBE"
    CONTINUOUS_STIRRED_TANK_REACTOR = "CONTINUOUS_STIRRED_TANK_REACTOR"
    PACKED_BED_REACTOR = "PACKED_BED_REACTOR"
    NMR_TUBE = "NMR_TUBE"
    PRESSURE_FLASK = "PRESSURE_FLASK"
    PRESSURE_REACTOR = "PRESSURE_REACTOR"
    ELECTROCHEMICAL_CELL = "ELECTROCHEMICAL_CELL"


class VesselMaterial(Enum):
    """ORD VesselMaterial.VesselMaterialType."""

    CUSTOM = "CUSTOM"
    GLASS = "GLASS"
    POLYPROPYLENE = "POLYPROPYLENE"
    PLASTIC = "PLASTIC"
    METAL = "METAL"
    QUARTZ = "QUARTZ"


class VesselPreparation(Enum):
    """ORD VesselPreparation.VesselPreparationType: how the vessel was prepared."""

    CUSTOM = "CUSTOM"
    NONE = "NONE"
    OVEN_DRIED = "OVEN_DRIED"
    FLAME_DRIED = "FLAME_DRIED"
    EVACUATED_BACKFILLED = "EVACUATED_BACKFILLED"
    PURGED = "PURGED"


class VesselAttachment(Enum):
    """ORD VesselAttachment.VesselAttachmentType: what is fitted to the vessel."""

    CUSTOM = "CUSTOM"
    SEPTUM = "SEPTUM"
    CAP = "CAP"
    MAT = "MAT"
    REFLUX_CONDENSER = "REFLUX_CONDENSER"
    VENT_NEEDLE = "VENT_NEEDLE"
    DEAN_STARK = "DEAN_STARK"
    VACUUM_TUBE = "VACUUM_TUBE"
    ADDITION_FUNNEL = "ADDITION_FUNNEL"
    DRYING_TUBE = "DRYING_TUBE"
    ALUMINUM_FOIL = "ALUMINUM_FOIL"
    THERMOCOUPLE = "THERMOCOUPLE"
    BALLOON = "BALLOON"
    GAS_ADAPTER = "GAS_ADAPTER"
    PRESSURE_REGULATOR = "PRESSURE_REGULATOR"
    RELEASE_VALVE = "RELEASE_VALVE"


class EnvironmentType(Enum):
    """ORD ReactionSetup.ReactionEnvironment: where the reaction is set up."""

    CUSTOM = "CUSTOM"
    FUME_HOOD = "FUME_HOOD"
    BENCH_TOP = "BENCH_TOP"
    GLOVE_BOX = "GLOVE_BOX"
    GLOVE_BAG = "GLOVE_BAG"


class TemperatureControl(Enum):
    """ORD TemperatureConditions.TemperatureControl: how the temperature is held."""

    CUSTOM = "CUSTOM"
    AMBIENT = "AMBIENT"
    OIL_BATH = "OIL_BATH"
    WATER_BATH = "WATER_BATH"
    SAND_BATH = "SAND_BATH"
    ICE_BATH = "ICE_BATH"
    DRY_ALUMINUM_PLATE = "DRY_ALUMINUM_PLATE"
    MICROWAVE = "MICROWAVE"
    DRY_ICE_BATH = "DRY_ICE_BATH"
    AIR_FAN = "AIR_FAN"
    LIQUID_NITROGEN = "LIQUID_NITROGEN"


class PressureControl(Enum):
    """ORD PressureConditions.PressureControl."""

    CUSTOM = "CUSTOM"
    AMBIENT = "AMBIENT"
    SLIGHT_POSITIVE = "SLIGHT_POSITIVE"
    SEALED = "SEALED"
    PRESSURIZED = "PRESSURIZED"


class Atmosphere(Enum):
    """ORD PressureConditions.Atmosphere: the gas above the reaction."""

    CUSTOM = "CUSTOM"
    AIR = "AIR"
    NITROGEN = "NITROGEN"
    ARGON = "ARGON"
    OXYGEN = "OXYGEN"
    HYDROGEN = "HYDROGEN"
    CARBON_MONOXIDE = "CARBON_MONOXIDE"
    CARBON_DIOXIDE = "CARBON_DIOXIDE"
    METHANE = "METHANE"
    AMMONIA = "AMMONIA"
    OZONE = "OZONE"
    ETHYLENE = "ETHYLENE"
    ACETYLENE = "ACETYLENE"


class StirringMethod(Enum):
    """ORD StirringConditions.StirringMethod."""

    CUSTOM = "CUSTOM"
    NONE = "NONE"
    STIR_BAR = "STIR_BAR"
    OVERHEAD_MIXER = "OVERHEAD_MIXER"
    AGITATION = "AGITATION"
    BALL_MILLING = "BALL_MILLING"
    SONICATION = "SONICATION"


class StirringRate(Enum):
    """ORD StirringConditions.StirringRate.StirringRateType, used when rpm is unknown."""

    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"


class AdditionSpeed(Enum):
    """ORD ReactionInput.AdditionSpeed.AdditionSpeedType."""

    ALL_AT_ONCE = "ALL_AT_ONCE"
    FAST = "FAST"
    SLOW = "SLOW"
    DROPWISE = "DROPWISE"
    CONTINUOUS = "CONTINUOUS"
    PORTIONWISE = "PORTIONWISE"


class AdditionDevice(Enum):
    """ORD ReactionInput.AdditionDevice.AdditionDeviceType."""

    CUSTOM = "CUSTOM"
    NONE = "NONE"
    SYRINGE = "SYRINGE"
    CANNULA = "CANNULA"
    ADDITION_FUNNEL = "ADDITION_FUNNEL"
    PIPETTE = "PIPETTE"
    POSITIVE_DISPLACEMENT_PIPETTE = "POSITIVE_DISPLACEMENT_PIPETTE"
    PISTON_PUMP = "PISTON_PUMP"
    SYRINGE_PUMP = "SYRINGE_PUMP"
    PERISTALTIC_PUMP = "PERISTALTIC_PUMP"


class AdditionForm(Enum):
    """
    How a material is added. Not an ORD enum: in ORD a solid or liquid is a
    ReactionInput component's texture, and a solution is a ReactionInput holding
    the solute and its solvent as two components.
    """

    SOLID = "SOLID"
    LIQUID = "LIQUID"
    SOLUTION = "SOLUTION"
    GAS = "GAS"


class MethodComponentRole(Enum):
    """Which part of the reaction table an added component comes from (ORD reaction_role)."""

    REACTANT = "REACTANT"
    REAGENT = "REAGENT"
    SOLVENT = "SOLVENT"
