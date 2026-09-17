from enum import Enum

from sources.extensions import db
from sources.models.section import Section


class WorkupStepType(Enum):
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


class WorkupSection(Section):
    """
    Workup section contains more information than default section
    """

    __tablename__ = "WorkupSection"

    id = db.Column(
        db.Integer,
        db.ForeignKey("Section.id", ondelete="CASCADE"),
        primary_key=True,
    )

    step_type = db.Column(Enum(WorkupStepType), nullable=False)

    details = db.Column(db.Text)

    __mapper_args__ = {"polymorphic_identity": "workup"}
