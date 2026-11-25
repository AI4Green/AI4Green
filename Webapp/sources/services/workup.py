"""
Reaction workups and purifications, built from typed steps, and the workbook
templates they can be saved as.

A reaction has an ordered list of WorkupSteps per category (WORKUP, PURIFICATION).
Saving them as a template keeps only the step types and their order; starting from a
template appends empty copies of those steps. Every change is sent to the audit log
through services.reaction_editing_history.
"""

from decimal import Decimal
from typing import Any, Dict, List, Optional

from flask_login import current_user
from sources import models, services
from sources.extensions import db
from sources.services.workup_steps import (
    WorkupValidationError,
    material_summary,
    validate_materials,
    validate_notes,
    validate_step_data,
    validate_temperature,
)


class WorkupNotFoundError(Exception):
    """Raised when a step or template doesn't belong to the reaction or workbook."""


def _record(
    reaction: models.Reaction, field_name: str, change_details: Dict[str, Any]
) -> None:
    """
    Sends an audit message for a workup change made by the current user.

    Args:
        reaction: The reaction whose workup changed.
        field_name: What happened.
        change_details: The details to record.
    """
    services.reaction_editing_history.workup_changed(
        services.person.from_current_user_email(),
        reaction,
        field_name,
        change_details,
    )


def parse_category(value: Any) -> models.WorkupCategory:
    """
    Reads a workup category from a request.

    Args:
        value: The submitted category.

    Returns:
        The category.

    Raises:
        WorkupValidationError: If it isn't a known category.
    """
    try:
        return models.WorkupCategory(value)
    except ValueError:
        raise WorkupValidationError("Unknown workup category.")


def parse_step_type(value: Any) -> models.WorkupStepType:
    """
    Reads a workup step type from a request.

    Args:
        value: The submitted step type.

    Returns:
        The step type.

    Raises:
        WorkupValidationError: If it isn't a known step type.
    """
    try:
        return models.WorkupStepType(value)
    except ValueError:
        raise WorkupValidationError("Unknown step type.")


# --- Steps --------------------------------------------------------------------


def list_steps(
    reaction: models.Reaction, category: models.WorkupCategory
) -> List[models.WorkupStep]:
    """
    Gets a reaction's steps for one category, in order.

    Args:
        reaction: The reaction.
        category: Workup or purification.

    Returns:
        The steps.
    """
    return (
        models.WorkupStep.query.filter_by(reaction_id=reaction.id, category=category)
        .order_by(models.WorkupStep.sort_order, models.WorkupStep.id)
        .all()
    )


def get_step(reaction: models.Reaction, step_id: int) -> models.WorkupStep:
    """
    Gets one of a reaction's steps.

    Args:
        reaction: The reaction.
        step_id: The step id.

    Returns:
        The step.

    Raises:
        WorkupNotFoundError: If the step doesn't exist or belongs to another reaction.
    """
    step = db.session.get(models.WorkupStep, step_id)
    if step is None or step.reaction_id != reaction.id:
        raise WorkupNotFoundError()
    return step


def _next_sort_order(reaction: models.Reaction, category: models.WorkupCategory) -> int:
    """
    Gets the position after a category's last step.

    Args:
        reaction: The reaction.
        category: Workup or purification.

    Returns:
        The next sort order.
    """
    highest = (
        db.session.query(db.func.max(models.WorkupStep.sort_order))
        .filter_by(reaction_id=reaction.id, category=category)
        .scalar()
    )
    return 1 if highest is None else highest + 1


def add_step(
    reaction: models.Reaction,
    category: models.WorkupCategory,
    step_type: models.WorkupStepType,
) -> models.WorkupStep:
    """
    Adds an empty step to the end of a reaction's workup or purification.

    Args:
        reaction: The reaction.
        category: Workup or purification.
        step_type: The kind of step.

    Returns:
        The new step.
    """
    step = models.WorkupStep(
        reaction_id=reaction.id,
        category=category,
        step_type=step_type,
        sort_order=_next_sort_order(reaction, category),
        data={},
        notes="",
    )
    db.session.add(step)
    db.session.commit()
    _record(
        reaction,
        "Workup Step Added",
        {"category": category.value, "step_type": step_type.value},
    )
    return step


def _audit_values(
    data: Dict[str, Any],
    notes: str,
    temperature: Optional[Decimal],
    materials: List[models.WorkupStepMaterial],
) -> Dict[str, Any]:
    """
    Flattens a step's values for get_differences, which doesn't handle nested dicts
    correctly. Each material slot becomes one key listing its materials.

    Args:
        data: The step's settings.
        notes: The step's notes.
        temperature: The step's temperature.
        materials: The step's materials.

    Returns:
        Flat key -> JSON-safe value.
    """
    values: Dict[str, Any] = {f"data.{key}": value for key, value in data.items()}
    values["notes"] = notes
    values["temperature"] = float(temperature) if temperature is not None else None
    for material in materials:
        values.setdefault(f"materials.{material.slot}", []).append(
            material_summary(material)
        )
    return values


def update_step(
    step: models.WorkupStep,
    data: Any,
    notes: Any,
    temperature: Any = None,
    materials: Any = None,
) -> models.WorkupStep:
    """
    Saves a step's settings, notes, temperature and materials, replacing what was
    there. Sends an audit message only if something changed.

    Args:
        step: The step.
        data: The submitted settings (validated against the step type).
        notes: The submitted notes.
        temperature: The submitted temperature in °C, or None.
        materials: The submitted materials (see workup_steps.validate_materials).

    Returns:
        The updated step.

    Raises:
        WorkupValidationError: If a value is invalid.
    """
    new_data = validate_step_data(step.step_type, data)
    new_notes = validate_notes(notes)
    new_temperature = validate_temperature(step.step_type, temperature)
    new_materials = validate_materials(
        step.step_type, step.reaction.workbook, materials
    )

    old = _audit_values(
        step.data or {}, step.notes or "", step.temperature, list(step.materials)
    )
    new = _audit_values(new_data, new_notes, new_temperature, new_materials)
    differences = services.reaction_editing_history.get_differences(old, new)
    if not differences:
        return step

    step.data = new_data
    step.notes = new_notes
    step.temperature = new_temperature
    step.materials = new_materials
    db.session.commit()
    _record(
        step.reaction,
        "Workup Step Edited",
        {
            "category": step.category.value,
            "step_type": step.step_type.value,
            "step_id": step.id,
            "changes": differences,
        },
    )
    return step


def reorder_steps(
    reaction: models.Reaction, category: models.WorkupCategory, step_ids: Any
) -> List[models.WorkupStep]:
    """
    Puts a category's steps in the given order.

    Args:
        reaction: The reaction.
        category: Workup or purification.
        step_ids: Every step id of that category, in the new order.

    Returns:
        The steps in their new order.

    Raises:
        WorkupValidationError: If the ids aren't exactly the category's steps.
    """
    steps = list_steps(reaction, category)
    by_id = {step.id: step for step in steps}
    if (
        not isinstance(step_ids, list)
        or len(step_ids) != len(by_id)
        or set(step_ids) != set(by_id)
    ):
        raise WorkupValidationError("Send every step of this list, in the new order.")

    for position, step_id in enumerate(step_ids, start=1):
        by_id[step_id].sort_order = position
    db.session.commit()
    _record(
        reaction,
        "Workup Steps Reordered",
        {
            "category": category.value,
            "order": [by_id[step_id].step_type.value for step_id in step_ids],
        },
    )
    return [by_id[step_id] for step_id in step_ids]


def delete_step(step: models.WorkupStep) -> None:
    """
    Deletes a step and closes the gap it leaves in the order.

    Args:
        step: The step.
    """
    reaction = step.reaction
    category = step.category
    details = {
        "category": category.value,
        "step_type": step.step_type.value,
        "data": step.data or {},
        "notes": step.notes or "",
        "temperature": (
            float(step.temperature) if step.temperature is not None else None
        ),
        "materials": [material_summary(material) for material in step.materials],
    }
    db.session.delete(step)
    db.session.flush()
    for position, remaining in enumerate(list_steps(reaction, category), start=1):
        remaining.sort_order = position
    db.session.commit()
    _record(reaction, "Workup Step Deleted", details)


# --- Templates ----------------------------------------------------------------


def list_templates(
    workbook: models.WorkBook, category: Optional[models.WorkupCategory] = None
) -> List[models.WorkupTemplate]:
    """
    Gets a workbook's workup templates, newest first.

    Args:
        workbook: The workbook.
        category: Only templates of this category, if given.

    Returns:
        The templates.
    """
    query = models.WorkupTemplate.query.filter_by(workbook_id=workbook.id)
    if category is not None:
        query = query.filter_by(category=category)
    return query.order_by(models.WorkupTemplate.time_of_creation.desc()).all()


def get_template(workbook: models.WorkBook, template_id: int) -> models.WorkupTemplate:
    """
    Gets one of a workbook's workup templates.

    Args:
        workbook: The workbook.
        template_id: The template id.

    Returns:
        The template.

    Raises:
        WorkupNotFoundError: If it doesn't exist or belongs to another workbook.
    """
    template = db.session.get(models.WorkupTemplate, template_id)
    if template is None or template.workbook_id != workbook.id:
        raise WorkupNotFoundError()
    return template


def save_as_template(
    reaction: models.Reaction,
    category: models.WorkupCategory,
    name: Any,
    description: Any,
) -> models.WorkupTemplate:
    """
    Saves a reaction's steps (types and order only) as a workbook template.

    Args:
        reaction: The reaction.
        category: Workup or purification.
        name: The template name.
        description: The template description.

    Returns:
        The new template.

    Raises:
        WorkupValidationError: If the name is missing or there are no steps.
    """
    name = name.strip() if isinstance(name, str) else ""
    if not name:
        raise WorkupValidationError("A template needs a name.")
    if len(name) > 200:
        raise WorkupValidationError("Template names must be 200 characters or fewer.")
    description = description.strip() if isinstance(description, str) else ""

    steps = list_steps(reaction, category)
    if not steps:
        raise WorkupValidationError("Add at least one step before saving a template.")

    template = models.WorkupTemplate(
        name=name,
        description=description,
        workbook_id=reaction.workbook.id,
        category=category,
        creator_id=current_user.id,
        steps=[
            models.WorkupTemplateStep(sort_order=position, step_type=step.step_type)
            for position, step in enumerate(steps, start=1)
        ],
    )
    db.session.add(template)
    db.session.commit()
    return template


def apply_template(
    reaction: models.Reaction, template: models.WorkupTemplate
) -> List[models.WorkupStep]:
    """
    Appends empty copies of a template's steps to the reaction's workup or
    purification (the template's category).

    Args:
        reaction: The reaction.
        template: The template.

    Returns:
        The reaction's steps for that category, in order.
    """
    next_position = _next_sort_order(reaction, template.category)
    for offset, template_step in enumerate(template.steps):
        db.session.add(
            models.WorkupStep(
                reaction_id=reaction.id,
                category=template.category,
                step_type=template_step.step_type,
                sort_order=next_position + offset,
                data={},
                notes="",
            )
        )
    db.session.commit()
    _record(
        reaction,
        "Workup Template Applied",
        {
            "category": template.category.value,
            "template_id": template.id,
            "template_name": template.name,
            "step_types": [step.step_type.value for step in template.steps],
        },
    )
    return list_steps(reaction, template.category)


def can_delete_template(template: models.WorkupTemplate) -> bool:
    """
    Whether the current user may delete a template (its creator).

    Args:
        template: The template.

    Returns:
        True for the template's creator.
    """
    return template.creator_id is not None and template.creator_id == current_user.id


def delete_template(template: models.WorkupTemplate) -> None:
    """
    Deletes a template. Reactions that used it keep their steps.

    Args:
        template: The template.
    """
    db.session.delete(template)
    db.session.commit()
