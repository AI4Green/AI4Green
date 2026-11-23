"""
What each workup step type records, and how its values are validated.

Each step type (StepDefinition) records up to three kinds of value:
- material slots (MaterialSlot): what was used, e.g. a wash's organic and aqueous
  phases. Each filled slot becomes a WorkupStepMaterial row linked to a Compound,
  NovelCompound or CHEM21 Solvent, so amounts can be queried for green chemistry
  metrics.
- a temperature in °C, saved in WorkupStep.temperature.
- other settings (StepField), saved in WorkupStep.data under fixed keys.

The definitions are also sent to the front end, which uses the labels, units and
choices to draw each step.
"""

import math
from dataclasses import dataclass, field
from decimal import Decimal
from typing import Any, Dict, List, Optional, Tuple

from flask_login import current_user
from sources import models, services
from sources.extensions import db
from sources.models import WorkupMaterialFate, WorkupMaterialRole, WorkupStepType

TEXT_MAX_LENGTH = 500
NOTES_MAX_LENGTH = 5000
MATERIAL_NAME_MAX_LENGTH = 200
MAX_MATERIALS_PER_STEP = 30
# the data key that repeated material slots are multiplied by
REPEATS_KEY = "repeats"
LOWEST_TEMPERATURE = -273.15
HIGHEST_TEMPERATURE = 1000
# largest amount Numeric(12, 4) can hold
MAX_QUANTITY = Decimal("99999999.9999")

TEXT = "text"
NUMBER = "number"
CHOICE = "choice"
# what the percentages of a mixture (e.g. an eluent) must add up to
PROPORTION_TOTAL = Decimal("100")

# what a filtration filters through; the front end only asks for the mass of
# filter aid (the "filter_aid" slot) for a Celite pad or silica plug
FILTER_MEDIA = (
    "Filter paper",
    "Sintered funnel",
    "Cotton wool",
    "Glass wool",
    "Celite pad",
    "Silica plug",
    "Other",
)

SOLVENT_PICKER = "solvent"
COMPOUND_PICKER = "compound"

# what the material added in an addition step is for
ADDITION_ROLES = (
    "Reagent",
    "Quench",
    "Internal standard",
    "Solvent",
    "Acid",
    "Base",
    "Antisolvent",
    "Catalyst",
    "Drying agent",
    "Seed crystal",
    "Other",
)

VOLUME_UNITS = ("mL", "L", "µL")
MASS_UNITS = ("g", "mg", "kg")
ANY_UNITS = VOLUME_UNITS + MASS_UNITS


class WorkupValidationError(Exception):
    """Raised when a workup step's values are invalid. The message is safe to show."""


@dataclass(frozen=True)
class StepField:
    """
    One setting a workup step records in WorkupStep.data.

    Attributes:
        key: Key in WorkupStep.data.
        label: Label shown to the user.
        kind: TEXT, NUMBER or CHOICE.
        unit: Unit shown next to a number, if any.
        choices: Allowed values of a CHOICE field.
        maximum: Largest allowed value of a NUMBER field, if any.
    """

    key: str
    label: str
    kind: str = TEXT
    unit: Optional[str] = None
    choices: Tuple[str, ...] = field(default_factory=tuple)
    maximum: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        """Serialises the field for the front end."""
        return {
            "key": self.key,
            "label": self.label,
            "kind": self.kind,
            "unit": self.unit,
            "choices": list(self.choices),
            "maximum": self.maximum,
        }


@dataclass(frozen=True)
class MaterialSlot:
    """
    A place in a step for a material, e.g. the wash's aqueous phase.

    Attributes:
        slot: Key saved in WorkupStepMaterial.slot.
        label: Label shown to the user.
        role: What the material does in the step, saved with it.
        picker: SOLVENT_PICKER (the coloured CHEM21 solvent list) or COMPOUND_PICKER
            (the compound database lookup).
        units: Allowed units; the first is the default.
        multiple: Whether the slot takes more than one material (e.g. eluent mixtures).
        repeated: Whether the amount is used once per `repeats` (e.g. each wash).
        proportions: Whether each material has a percentage of the mixture (e.g.
            90% hexane, 10% EtOAc in an eluent). They must add up to 100.
        default_fate: What happens to the material unless the user says otherwise.
    """

    slot: str
    label: str
    role: WorkupMaterialRole
    picker: str
    units: Tuple[str, ...]
    multiple: bool = False
    repeated: bool = False
    proportions: bool = False
    default_fate: WorkupMaterialFate = WorkupMaterialFate.WASTE

    def to_dict(self) -> Dict[str, Any]:
        """Serialises the slot for the front end."""
        return {
            "slot": self.slot,
            "label": self.label,
            "role": self.role.value,
            "picker": self.picker,
            "units": list(self.units),
            "multiple": self.multiple,
            "repeated": self.repeated,
            "proportions": self.proportions,
            "defaultFate": self.default_fate.value,
        }


@dataclass(frozen=True)
class StepDefinition:
    """
    Everything one step type records.

    Attributes:
        fields: Settings saved in WorkupStep.data.
        materials: Material slots, saved as WorkupStepMaterial rows.
        temperature_label: Label of the step's temperature, or None if it has none.
    """

    fields: Tuple[StepField, ...] = ()
    materials: Tuple[MaterialSlot, ...] = ()
    temperature_label: Optional[str] = None

    def slot(self, key: str) -> Optional[MaterialSlot]:
        """
        Finds one of the step type's material slots.

        Args:
            key: The slot key.

        Returns:
            The slot, or None if the step type doesn't have it.
        """
        return next((slot for slot in self.materials if slot.slot == key), None)

    def to_dict(self) -> Dict[str, Any]:
        """Serialises the definition for the front end."""
        return {
            "fields": [step_field.to_dict() for step_field in self.fields],
            "materials": [slot.to_dict() for slot in self.materials],
            "temperature": self.temperature_label,
        }


def _text(key: str, label: str) -> StepField:
    """Builds a text field."""
    return StepField(key=key, label=label)


def _number(
    key: str, label: str, unit: Optional[str] = None, maximum: Optional[float] = None
) -> StepField:
    """Builds a number field."""
    return StepField(key=key, label=label, kind=NUMBER, unit=unit, maximum=maximum)


def _choice(key: str, label: str, choices: Tuple[str, ...]) -> StepField:
    """Builds a choice field."""
    return StepField(key=key, label=label, kind=CHOICE, choices=choices)


def _repeats(label: str) -> StepField:
    """Builds the field that repeated material slots are multiplied by."""
    return _number(REPEATS_KEY, label)


def _solvent(
    slot: str,
    label: str,
    role: WorkupMaterialRole = WorkupMaterialRole.SOLVENT,
    **options: Any,
) -> MaterialSlot:
    """Builds a slot picked from the CHEM21 solvent list, measured by volume."""
    return MaterialSlot(
        slot=slot,
        label=label,
        role=role,
        picker=SOLVENT_PICKER,
        units=VOLUME_UNITS,
        **options,
    )


def _compound(
    slot: str,
    label: str,
    role: WorkupMaterialRole,
    units: Tuple[str, ...] = MASS_UNITS,
    **options: Any,
) -> MaterialSlot:
    """Builds a slot picked from the compound database."""
    return MaterialSlot(
        slot=slot,
        label=label,
        role=role,
        picker=COMPOUND_PICKER,
        units=units,
        **options,
    )


Role = WorkupMaterialRole

STEP_DEFINITIONS: Dict[WorkupStepType, StepDefinition] = {
    WorkupStepType.CUSTOM: StepDefinition(
        fields=(_text("title", "Step name"),),
        materials=(
            _compound(
                "materials", "Materials used", Role.REAGENT, ANY_UNITS, multiple=True
            ),
        ),
    ),
    WorkupStepType.ADDITION: StepDefinition(
        fields=(
            _text("rate", "Addition rate"),
            _choice("role", "Role", ADDITION_ROLES),
        ),
        materials=(_compound("material", "Material added", Role.REAGENT, ANY_UNITS),),
        temperature_label="Temperature",
    ),
    WorkupStepType.ALIQUOT: StepDefinition(
        fields=(
            _number("amount", "Amount"),
            _choice("amount_unit", "Unit", ("mL", "µL", "%")),
            _text("purpose", "Destination or purpose"),
        ),
    ),
    WorkupStepType.TEMPERATURE: StepDefinition(
        fields=(_text("duration", "Duration"),),
        temperature_label="Target temperature",
    ),
    WorkupStepType.CONCENTRATION: StepDefinition(
        fields=(
            _choice(
                "method", "Method", ("Rotary evaporator", "Nitrogen stream", "Vacuum")
            ),
            _number("pressure", "Pressure", "mbar"),
        ),
        temperature_label="Bath temperature",
    ),
    WorkupStepType.EXTRACTION: StepDefinition(
        fields=(
            _repeats("Number of extractions"),
            _choice("phase_kept", "Phase kept", ("Organic", "Aqueous")),
        ),
        materials=(_solvent("solvent", "Extraction solvent", repeated=True),),
    ),
    WorkupStepType.FILTRATION: StepDefinition(
        fields=(
            _choice("method", "Type", ("Vacuum", "Gravity")),
            _choice("filter_medium", "Filter medium", FILTER_MEDIA),
            _choice("phase_kept", "Kept", ("Filtrate", "Solid")),
        ),
        materials=(
            _compound("filter_aid", "Celite or silica used", Role.REAGENT),
            _solvent("rinse_solvent", "Rinse solvent"),
        ),
    ),
    WorkupStepType.WASH: StepDefinition(
        fields=(_repeats("Number of washes"),),
        materials=(
            _solvent(
                "organic_phase",
                "Organic phase",
                Role.ORGANIC_PHASE,
                default_fate=WorkupMaterialFate.KEPT,
            ),
            _compound(
                "aqueous_phase",
                "Aqueous phase",
                Role.AQUEOUS_PHASE,
                VOLUME_UNITS,
                repeated=True,
            ),
        ),
    ),
    WorkupStepType.DRY_IN_VACUUM: StepDefinition(
        fields=(
            _number("pressure", "Pressure", "mbar"),
            _text("duration", "Duration"),
        ),
        temperature_label="Temperature",
    ),
    WorkupStepType.DRY_WITH_MATERIAL: StepDefinition(
        materials=(_compound("drying_agent", "Drying agent", Role.DRYING_AGENT),),
    ),
    WorkupStepType.SCAVENGING: StepDefinition(
        fields=(_text("duration", "Duration"),),
        materials=(_compound("scavenger", "Scavenger", Role.REAGENT),),
    ),
    WorkupStepType.WAIT: StepDefinition(fields=(_text("duration", "Duration"),)),
    WorkupStepType.STIRRING: StepDefinition(
        fields=(
            _text("duration", "Duration"),
            _number("rate", "Stirring rate", "rpm"),
        ),
    ),
    WorkupStepType.PH_ADJUST: StepDefinition(
        fields=(_number("target_ph", "Target pH", maximum=14),),
        materials=(_compound("reagent", "Reagent used", Role.REAGENT, ANY_UNITS),),
    ),
    WorkupStepType.DISSOLUTION: StepDefinition(
        materials=(_solvent("solvent", "Solvent"),),
    ),
    WorkupStepType.FLASH_CHROMATOGRAPHY: StepDefinition(
        fields=(_text("gradient", "Gradient"),),
        materials=(
            _compound("stationary_phase", "Stationary phase", Role.STATIONARY_PHASE),
            _solvent("eluent", "Eluent", Role.ELUENT, multiple=True, proportions=True),
        ),
    ),
    WorkupStepType.OTHER_CHROMATOGRAPHY: StepDefinition(
        fields=(
            _choice(
                "method",
                "Method",
                ("Preparative TLC", "Preparative HPLC", "Column", "Other"),
            ),
        ),
        materials=(
            _compound("stationary_phase", "Stationary phase", Role.STATIONARY_PHASE),
            _solvent("eluent", "Eluent", Role.ELUENT, multiple=True, proportions=True),
        ),
    ),
    WorkupStepType.DISTILLATION: StepDefinition(
        fields=(
            _choice(
                "method", "Method", ("Simple", "Fractional", "Vacuum", "Kugelrohr")
            ),
            _number("pressure", "Pressure", "mbar"),
        ),
        temperature_label="Temperature",
    ),
    WorkupStepType.CRYSTALLISATION: StepDefinition(
        fields=(_text("duration", "Duration"),),
        materials=(_solvent("solvent", "Solvent", multiple=True),),
        temperature_label="Temperature",
    ),
}


def definitions_as_dict() -> Dict[str, Dict[str, Any]]:
    """
    Serialises the step definitions for the front end.

    Returns:
        Step type value -> `{fields, materials, temperature}`.
    """
    return {
        step_type.value: definition.to_dict()
        for step_type, definition in STEP_DEFINITIONS.items()
    }


# --- Settings (WorkupStep.data) -----------------------------------------------


def _finite_number(label: str, value: Any) -> float:
    """
    Reads a finite number from a submitted value.

    Args:
        label: The value's label, for the error message.
        value: The submitted value.

    Returns:
        The number.

    Raises:
        WorkupValidationError: If it isn't a finite number.
    """
    if isinstance(value, bool):
        raise WorkupValidationError(f"{label} must be a number.")
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise WorkupValidationError(f"{label} must be a number.")
    if not math.isfinite(number):
        raise WorkupValidationError(f"{label} must be a number.")
    return number


def _clean_number(step_field: StepField, value: Any) -> float:
    """
    Validates a number value.

    Args:
        step_field: The field definition.
        value: The submitted value.

    Returns:
        The value as a number.
    """
    number = _finite_number(step_field.label, value)
    if number < 0:
        raise WorkupValidationError(f"{step_field.label} must be zero or more.")
    if step_field.maximum is not None and number > step_field.maximum:
        raise WorkupValidationError(
            f"{step_field.label} must be {step_field.maximum:g} or less."
        )
    return int(number) if number.is_integer() else number


def _clean_text(step_field: StepField, value: Any) -> str:
    """
    Validates a text value.

    Args:
        step_field: The field definition.
        value: The submitted value.

    Returns:
        The trimmed text.
    """
    if not isinstance(value, str):
        raise WorkupValidationError(f"{step_field.label} must be text.")
    text = value.strip()
    if len(text) > TEXT_MAX_LENGTH:
        raise WorkupValidationError(
            f"{step_field.label} must be {TEXT_MAX_LENGTH} characters or fewer."
        )
    return text


def _clean_choice(step_field: StepField, value: Any) -> str:
    """
    Validates a choice value.

    Args:
        step_field: The field definition.
        value: The submitted value.

    Returns:
        The choice.
    """
    if value not in step_field.choices:
        raise WorkupValidationError(
            f"{step_field.label} must be one of: {', '.join(step_field.choices)}."
        )
    return value


def validate_step_data(step_type: WorkupStepType, data: Any) -> Dict[str, Any]:
    """
    Cleans a step's submitted settings against its definition.

    Unknown keys are dropped and empty values are left out, so the saved data only
    holds what the user entered.

    Args:
        step_type: The step's type.
        data: The submitted settings, a JSON object.

    Returns:
        The cleaned settings.

    Raises:
        WorkupValidationError: If the data isn't an object or a value is invalid.
    """
    if not isinstance(data, dict):
        raise WorkupValidationError("Step values must be a JSON object.")

    cleaners = {NUMBER: _clean_number, TEXT: _clean_text, CHOICE: _clean_choice}
    cleaned: Dict[str, Any] = {}
    for step_field in STEP_DEFINITIONS[step_type].fields:
        value = data.get(step_field.key)
        if value is None or (isinstance(value, str) and not value.strip()):
            continue
        cleaned[step_field.key] = cleaners[step_field.kind](step_field, value)
    return cleaned


def validate_temperature(step_type: WorkupStepType, value: Any) -> Optional[Decimal]:
    """
    Validates a step's temperature.

    Args:
        step_type: The step's type.
        value: The submitted temperature in °C, or None/"" for none.

    Returns:
        The temperature rounded to 0.01 °C, or None.

    Raises:
        WorkupValidationError: If the step type has no temperature or it's out of range.
    """
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    label = STEP_DEFINITIONS[step_type].temperature_label
    if label is None:
        raise WorkupValidationError("This step doesn't record a temperature.")
    number = _finite_number(label, value)
    if not LOWEST_TEMPERATURE <= number <= HIGHEST_TEMPERATURE:
        raise WorkupValidationError(
            f"{label} must be between {LOWEST_TEMPERATURE:g} and"
            f" {HIGHEST_TEMPERATURE:g} °C."
        )
    return Decimal(str(number)).quantize(Decimal("0.01"))


def validate_notes(notes: Any) -> str:
    """
    Validates a step's free-text notes.

    Args:
        notes: The submitted notes.

    Returns:
        The trimmed notes.

    Raises:
        WorkupValidationError: If the notes aren't text or are too long.
    """
    if notes is None:
        return ""
    if not isinstance(notes, str):
        raise WorkupValidationError("Notes must be text.")
    if len(notes) > NOTES_MAX_LENGTH:
        raise WorkupValidationError(
            f"Notes must be {NOTES_MAX_LENGTH} characters or fewer."
        )
    return notes.strip()


# --- Materials (WorkupStepMaterial) -------------------------------------------


def _optional_id(label: str, value: Any) -> Optional[int]:
    """
    Reads an optional positive integer id.

    Args:
        label: What the id is for, for the error message.
        value: The submitted id.

    Returns:
        The id, or None.

    Raises:
        WorkupValidationError: If it isn't a positive integer.
    """
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise WorkupValidationError(f"{label} isn't a valid id.")
    return value


def _accessible_workbook_ids(workbook: models.WorkBook) -> List[int]:
    """
    The workbooks whose novel compounds the current user may link: theirs in the
    reaction's workgroup, the same scope as the compound lookup.

    Args:
        workbook: The reaction's workbook.

    Returns:
        The workbook ids.
    """
    return [
        wb.id
        for wb in services.workbook.get_workbooks_from_user_group_combination(
            workbook.WorkGroup.name
        )
    ]


def _accessible_solvent(
    workbook: models.WorkBook, name: str
) -> Optional[models.Solvent]:
    """
    Finds a solvent on the workbook's solvent list (CHEM21 solvents plus the
    workbook's own novel solvents).

    Args:
        workbook: The reaction's workbook.
        name: The solvent name.

    Returns:
        The solvent, or None if it isn't on the list.
    """
    solvent = db.session.get(models.Solvent, name)
    if solvent is None:
        return None
    if solvent.novel_compound and solvent.novel_compound[0].workbook != workbook.id:
        return None
    return solvent


def _linked_name(
    workbook: models.WorkBook,
    compound_id: Optional[int],
    novel_compound_id: Optional[int],
    solvent_name: Optional[str],
) -> Optional[str]:
    """
    Checks a material's link points at something the user may use and returns its
    name.

    Args:
        workbook: The reaction's workbook.
        compound_id: A Compound id, or None.
        novel_compound_id: A NovelCompound id, or None.
        solvent_name: A Solvent name, or None.

    Returns:
        The linked record's name, or None if the material isn't linked.

    Raises:
        WorkupValidationError: If more than one link is given or the link isn't usable.
    """
    links = [x for x in (compound_id, novel_compound_id, solvent_name) if x is not None]
    if len(links) > 1:
        raise WorkupValidationError("A material can only be linked to one compound.")
    if compound_id is not None:
        compound = db.session.get(models.Compound, compound_id)
        if compound is None:
            raise WorkupValidationError("That compound couldn't be found.")
        return compound.name
    if novel_compound_id is not None:
        novel_compound = db.session.get(models.NovelCompound, novel_compound_id)
        if novel_compound is None or novel_compound.workbook not in (
            _accessible_workbook_ids(workbook)
        ):
            raise WorkupValidationError("That compound couldn't be found.")
        return novel_compound.name
    if solvent_name is not None:
        solvent = _accessible_solvent(workbook, solvent_name)
        if solvent is None:
            raise WorkupValidationError("That solvent couldn't be found.")
        return solvent.name
    return None


def _clean_quantity(label: str, value: Any) -> Optional[Decimal]:
    """
    Validates a material's amount.

    Args:
        label: The slot label, for the error message.
        value: The submitted amount, or None/"" for an unknown amount.

    Returns:
        The amount, or None.
    """
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    number = _finite_number(f"{label} amount", value)
    if number < 0:
        raise WorkupValidationError(f"{label} amount must be zero or more.")
    quantity = Decimal(str(number)).quantize(Decimal("0.0001"))
    if quantity > MAX_QUANTITY:
        raise WorkupValidationError(f"{label} amount is too large.")
    return quantity


def _clean_proportion(label: str, value: Any) -> Optional[Decimal]:
    """
    Validates a material's percentage of its mixture.

    Args:
        label: The slot label, for the error message.
        value: The submitted percentage, or None/"" for none.

    Returns:
        The percentage rounded to 0.01, or None.
    """
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    number = _finite_number(f"{label} percentage", value)
    if not 0 <= number <= PROPORTION_TOTAL:
        raise WorkupValidationError(
            f"{label} percentage must be between 0 and {PROPORTION_TOTAL}."
        )
    return Decimal(str(number)).quantize(Decimal("0.01"))


def _check_proportions(
    slot: MaterialSlot, rows: List[models.WorkupStepMaterial]
) -> None:
    """
    Checks that the percentages of a slot's materials are all given, or none are,
    and that they add up to 100.

    Args:
        slot: A slot with proportions.
        rows: The slot's cleaned materials.

    Raises:
        WorkupValidationError: If only some are given or they don't add up to 100.
    """
    given = [row.proportion for row in rows if row.proportion is not None]
    if not given:
        return
    if len(given) < len(rows):
        raise WorkupValidationError(
            f"Give a percentage for each {slot.label.lower()}, or for none."
        )
    total = sum(given, Decimal("0"))
    if total != PROPORTION_TOTAL:
        raise WorkupValidationError(
            f"{slot.label} percentages must add up to {PROPORTION_TOTAL}%"
            f" (they add up to {total.normalize():f}%)."
        )


def _clean_material(
    workbook: models.WorkBook,
    slot: MaterialSlot,
    material: Dict[str, Any],
    sort_order: int,
) -> Optional[models.WorkupStepMaterial]:
    """
    Validates one submitted material and builds its row.

    Args:
        workbook: The reaction's workbook.
        slot: The slot it fills.
        material: The submitted material.
        sort_order: Its position within the slot.

    Returns:
        The row, or None if the material is completely empty.
    """
    compound_id = _optional_id("Compound", material.get("compoundId"))
    novel_compound_id = _optional_id("Compound", material.get("novelCompoundId"))
    solvent_name = material.get("solventName") or None
    if solvent_name is not None and not isinstance(solvent_name, str):
        raise WorkupValidationError("That solvent couldn't be found.")

    name = material.get("name") or ""
    if not isinstance(name, str):
        raise WorkupValidationError(f"{slot.label} must be text.")
    name = name.strip()
    if len(name) > MATERIAL_NAME_MAX_LENGTH:
        raise WorkupValidationError(
            f"{slot.label} must be {MATERIAL_NAME_MAX_LENGTH} characters or fewer."
        )

    linked_name = _linked_name(workbook, compound_id, novel_compound_id, solvent_name)
    quantity = _clean_quantity(slot.label, material.get("quantity"))
    if linked_name is None and not name and quantity is None:
        return None

    unit = material.get("unit") or slot.units[0]
    if unit not in slot.units:
        raise WorkupValidationError(
            f"{slot.label} unit must be one of: {', '.join(slot.units)}."
        )
    try:
        fate = WorkupMaterialFate(material.get("fate") or slot.default_fate.value)
    except ValueError:
        raise WorkupValidationError("Unknown material fate.")
    proportion = (
        _clean_proportion(slot.label, material.get("proportion"))
        if slot.proportions
        else None
    )

    return models.WorkupStepMaterial(
        slot=slot.slot,
        sort_order=sort_order,
        role=slot.role,
        fate=fate,
        compound_id=compound_id,
        novel_compound_id=novel_compound_id,
        solvent_name=solvent_name,
        name=linked_name or name,
        quantity=quantity,
        unit=unit,
        proportion=proportion,
    )


def validate_materials(
    step_type: WorkupStepType, workbook: models.WorkBook, materials: Any
) -> List[models.WorkupStepMaterial]:
    """
    Cleans a step's submitted materials against its slots.

    Each material is `{slot, compoundId?, novelCompoundId?, solventName?, name,
    quantity, unit, fate, proportion?}`. A linked material takes its name from the
    linked record; completely empty materials are dropped. `proportion` (a
    percentage) is only kept for slots with proportions, where the slot's
    percentages must all be given (or none) and add up to 100.

    Args:
        step_type: The step's type.
        workbook: The reaction's workbook, which limits the novel compounds and
            solvents that can be linked.
        materials: The submitted materials, a list (None for none).

    Returns:
        New, unsaved WorkupStepMaterial rows.

    Raises:
        WorkupValidationError: If a material is invalid or doesn't fit its slot.
    """
    if materials is None:
        return []
    if not isinstance(materials, list):
        raise WorkupValidationError("Materials must be a list.")
    if len(materials) > MAX_MATERIALS_PER_STEP:
        raise WorkupValidationError(
            f"A step can have at most {MAX_MATERIALS_PER_STEP} materials."
        )

    definition = STEP_DEFINITIONS[step_type]
    rows: List[models.WorkupStepMaterial] = []
    count_by_slot: Dict[str, int] = {}
    for material in materials:
        if not isinstance(material, dict):
            raise WorkupValidationError("Each material must be a JSON object.")
        slot = definition.slot(material.get("slot"))
        if slot is None:
            raise WorkupValidationError("This step doesn't have that material.")
        position = count_by_slot.get(slot.slot, 0) + 1
        row = _clean_material(workbook, slot, material, position)
        if row is None:
            continue
        if position > 1 and not slot.multiple:
            raise WorkupValidationError(f"{slot.label} takes one material.")
        count_by_slot[slot.slot] = position
        rows.append(row)

    for slot in definition.materials:
        if slot.proportions:
            _check_proportions(slot, [row for row in rows if row.slot == slot.slot])
    return rows


def total_quantity(
    material: models.WorkupStepMaterial,
) -> Tuple[Optional[Decimal], str]:
    """
    The total amount of a material used in its step. Amounts in repeated slots
    (e.g. each wash) are multiplied by the step's number of repeats. Green chemistry
    metrics must use this rather than `quantity`.

    Args:
        material: The material.

    Returns:
        The total amount (None if unknown) and its unit.
    """
    if material.quantity is None:
        return None, material.unit

    step = material.step
    slot = STEP_DEFINITIONS[step.step_type].slot(material.slot)
    repeats = (step.data or {}).get(REPEATS_KEY) if slot and slot.repeated else None
    if not repeats:
        return material.quantity, material.unit
    return material.quantity * Decimal(str(repeats)), material.unit


def material_summary(material: models.WorkupStepMaterial) -> str:
    """
    A one-line description of a material for the audit log, e.g.
    "Testanol 20 mL (waste)" or "Hexane 90% 450 mL (waste)". Saving a step uses
    it to spot changes, so it includes everything the user can edit.

    Args:
        material: The material.

    Returns:
        The description.
    """
    share = (
        f" {material.proportion.normalize():f}%"
        if material.proportion is not None
        else ""
    )
    amount = (
        f" {material.quantity.normalize():f} {material.unit}"
        if material.quantity is not None
        else ""
    )
    name = material.name or "(unnamed)"
    return f"{name}{share}{amount} ({material.fate.value.lower()})"


# --- Lookups ------------------------------------------------------------------


def solvent_options(workbook: models.WorkBook) -> List[Dict[str, Any]]:
    """
    The workbook's solvent list (the same one as the reaction table), with each
    solvent's CHEM21 flag and the record a material picked from it should link to.

    Args:
        workbook: The reaction's workbook.

    Returns:
        `{name, flag, compoundId, novelCompoundId, solventName}` per solvent, by name.
        Exactly one of the three links is set.
    """
    options = []
    for solvent in services.solvent.get_workbook_list(workbook):
        compound = solvent.compound[0] if solvent.compound else None
        novel_compound = solvent.novel_compound[0] if solvent.novel_compound else None
        options.append(
            {
                "name": solvent.name,
                "flag": solvent.flag,
                "compoundId": compound.id if compound else None,
                "novelCompoundId": (
                    novel_compound.id if novel_compound and not compound else None
                ),
                "solventName": (
                    solvent.name
                    if compound is None and novel_compound is None
                    else None
                ),
            }
        )
    return sorted(options, key=lambda option: option["name"].lower())


def hazard_colours() -> Dict[str, str]:
    """
    The current user's CHEM21 colours (they can change them in their account), so
    the solvent list looks the same as on the reaction table.

    Returns:
        Colour name -> hex colour, e.g. {"Recommended": "#00ff00", ...}.
    """
    return dict(current_user.hazard_colors or {})
