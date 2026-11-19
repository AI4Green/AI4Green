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
