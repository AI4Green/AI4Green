from enum import Enum
from typing import Any, Dict, Optional

from sources.extensions import db
from sources.models.section import Section


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


class WorkupSection(Section):
    """
    One workup step. Adds the ORD step type and free-text details to a section.
    """

    __tablename__ = "WorkupSection"

    id = db.Column(
        db.Integer,
        db.ForeignKey("Section.id", ondelete="CASCADE"),
        primary_key=True,
    )

    step_type = db.Column(db.Enum(WorkupStepType), nullable=False)

    details = db.Column(db.Text)

    __mapper_args__ = {"polymorphic_identity": "workup_section"}

    def to_dict(self, template_instance_id: Optional[int] = None) -> Dict[str, Any]:
        """
        Serialises the step: the section fields plus step type and details.

        Args:
            template_instance_id: If given, include that instance's responses and added fields.

        Returns:
            Dictionary representation of the step.
        """
        data = super().to_dict(template_instance_id)
        data.update({"stepType": self.step_type.value, "details": self.details})
        return data
