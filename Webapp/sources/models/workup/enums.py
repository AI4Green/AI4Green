from enum import Enum


class WorkupCategory(Enum):
    """
    Which reaction card a workup step or template belongs to.

    ORD has a single ReactionWorkup type covering both; this only splits them for users.
    """

    WORKUP = "WORKUP"
    PURIFICATION = "PURIFICATION"


class WorkupStepType(Enum):
    """
    Workup step types, from the Open Reaction Database ReactionWorkup.WorkupType.
    ORD treats purification (chromatography, distillation, crystallisation) as workup too.
    """

    CUSTOM = "CUSTOM"
    ADDITION = "ADDITION"
    ALIQUOT = "ALIQUOT"
    TEMPERATURE = "TEMPERATURE"
    CONCENTRATION = "CONCENTRATION"
    EXTRACTION = "EXTRACTION"
    FILTRATION = "FILTRATION"
    WASH = "WASH"
    DRY_IN_VACUUM = "DRY_IN_VACUUM"
    DRY_WITH_MATERIAL = "DRY_WITH_MATERIAL"
    SCAVENGING = "SCAVENGING"
    WAIT = "WAIT"
    STIRRING = "STIRRING"
    PH_ADJUST = "PH_ADJUST"
    DISSOLUTION = "DISSOLUTION"
    FLASH_CHROMATOGRAPHY = "FLASH_CHROMATOGRAPHY"
    OTHER_CHROMATOGRAPHY = "OTHER_CHROMATOGRAPHY"
    DISTILLATION = "DISTILLATION"
    CRYSTALLISATION = "CRYSTALLISATION"


class WorkupMaterialRole(Enum):
    """What a material is doing in a workup step, for green chemistry metrics."""

    SOLVENT = "SOLVENT"
    ORGANIC_PHASE = "ORGANIC_PHASE"
    AQUEOUS_PHASE = "AQUEOUS_PHASE"
    REAGENT = "REAGENT"
    QUENCH = "QUENCH"
    DRYING_AGENT = "DRYING_AGENT"
    STATIONARY_PHASE = "STATIONARY_PHASE"
    ELUENT = "ELUENT"
    INTERNAL_STANDARD = "INTERNAL_STANDARD"


class WorkupMaterialFate(Enum):
    """
    What happens to a material after the step: discarded as waste, recovered for
    reuse, or kept (e.g. the organic phase carried on to the next step).
    """

    WASTE = "WASTE"
    RECOVERED = "RECOVERED"
    KEPT = "KEPT"
