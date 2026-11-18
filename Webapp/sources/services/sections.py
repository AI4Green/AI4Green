"""
Sections and their fields, for templates and for template instances.

A section belongs to a template, or to one template instance that added it (see
models.Section). Fields on a template section are shared by every instance; an
instance can add its own fields to a template section (Field.template_instance_id).
"""

from typing import Any, Dict, List, Optional

from sources import models
from sources.extensions import db


class SectionInUseError(Exception):
    """Raised when removing a section or field that already has saved values."""


def _template_needs_workup_sections(template: models.Template) -> bool:
    """
    Whether sections of this template must be workup steps.

    Args:
        template: The template that will own the section.

    Returns:
        True for workup templates.
    """
    return isinstance(template, models.WorkupTemplate)


def create_section(
    name: str,
    sort_order: Optional[int],
    template: Optional[models.Template] = None,
    template_instance: Optional[models.WorkupInstance] = None,
    step_type: Optional[str] = None,
    details: Optional[str] = None,
    commit: bool = True,
) -> models.Section:
    """
    Creates a section on a template, or on one template instance.

    Workup templates and workup instances always get a WorkupSection (step type
    CUSTOM if none is given). Other templates get a plain Section.

    Args:
        name: Section name.
        sort_order: Position of the section.
        template: The owning template. Give this or template_instance, not both.
        template_instance: The owning instance, for a section only that instance has.
        step_type: ORD step type (models.WorkupStepType value), for workup sections.
        details: Free-text step details, for workup sections.
        commit: Whether to commit the session.

    Returns:
        The new section.
    """
    if (template is None) == (template_instance is None):
        raise ValueError("A section needs exactly one owner")

    owner_template = template or template_instance.template
    owner_kwargs = (
        {"template_id": template.id}
        if template is not None
        else {"template_instance_id": template_instance.id}
    )

    if _template_needs_workup_sections(owner_template):
        section = models.WorkupSection(
            name=name,
            sort_order=sort_order,
            step_type=models.WorkupStepType(step_type or "CUSTOM"),
            details=details,
            **owner_kwargs,
        )
    else:
        section = models.Section(name=name, sort_order=sort_order, **owner_kwargs)

    db.session.add(section)
    if commit:
        db.session.commit()
    return section


def update_template_sections(
    template: models.Template, sections_data: List[Dict[str, Any]]
) -> List[models.Section]:
    """
    Renames and reorders a template's sections, and creates any new ones.

    Section ids that don't belong to this template are ignored, so one template's
    request can't change another template.

    Args:
        template: The template being edited.
        sections_data: Items with "name", "sortOrder" and, for existing sections, "id".
            Workup items can also have "stepType" and "details".

    Returns:
        The updated and created sections.
    """
    own_sections = {section.id: section for section in template.sections}
    updated = []

    for item in sections_data:
        section_id = item.get("id")
        if section_id is not None:
            section = own_sections.get(section_id)
            if section is None:
                continue
            section.name = item.get("name")
            section.sort_order = item.get("sortOrder")
            if isinstance(section, models.WorkupSection):
                if item.get("stepType"):
                    section.step_type = models.WorkupStepType(item["stepType"])
                if "details" in item:
                    section.details = item.get("details")
            updated.append(section)
        else:
            updated.append(
                create_section(
                    name=item.get("name"),
                    sort_order=item.get("sortOrder"),
                    template=template,
                    step_type=item.get("stepType"),
                    details=item.get("details"),
                    commit=False,
                )
            )

    db.session.commit()
    return updated


def editable_fields(
    section: models.Section,
    template_instance: Optional[models.WorkupInstance] = None,
) -> List[models.Field]:
    """
    The fields of a section that can be changed in this context.

    - Editing a template: the template's own fields.
    - An instance editing a section it added: all of that section's fields.
    - An instance editing a template section: only the fields it added there.

    Args:
        section: The section being edited.
        template_instance: The instance doing the editing, or None for the template.

    Returns:
        The editable fields.
    """
    if template_instance is None or section.template_instance_id is not None:
        return [field for field in section.fields if field.template_instance_id is None]
    return [
        field
        for field in section.fields
        if field.template_instance_id == template_instance.id
    ]


def _field_has_responses(field_ids: List[int]) -> bool:
    """
    Whether any of these fields has a saved response.

    Args:
        field_ids: Field ids to check.

    Returns:
        True if at least one response exists.
    """
    if not field_ids:
        return False
    return (
        db.session.query(models.FieldResponse.id)
        .filter(models.FieldResponse.field_id.in_(field_ids))
        .first()
        is not None
    )


def save_fields(
    section: models.Section,
    fields_data: List[Dict[str, Any]],
    template_instance: Optional[models.WorkupInstance] = None,
) -> None:
    """
    Replaces the editable fields of a section with the submitted list.

    Fields missing from the list are deleted, fields with an id are updated and
    fields without one are created. Only fields returned by editable_fields can be
    touched, so an instance can never change or delete a template field.

    Args:
        section: The section being edited.
        fields_data: Items with "name", "sortOrder", "inputType" (InputType id),
            "mandatory" and, for existing fields, "id".
        template_instance: The instance doing the editing, or None for the template.

    Raises:
        SectionInUseError: When a field to delete has saved values in any instance.
    """
    existing = {
        field.id: field for field in editable_fields(section, template_instance)
    }
    incoming_ids = {item.get("id") for item in fields_data if item.get("id")}

    ids_to_delete = [field_id for field_id in existing if field_id not in incoming_ids]
    # saved values are kept for the audit trail, so a field with any can't go
    if _field_has_responses(ids_to_delete):
        raise SectionInUseError("This field has saved values, so it can't be deleted.")
    for field_id in ids_to_delete:
        db.session.delete(existing[field_id])

    # fields an instance adds to a template section are only shown to that instance
    new_field_instance_id = (
        template_instance.id
        if template_instance is not None and section.template_instance_id is None
        else None
    )

    for item in fields_data:
        values = {
            "name": item.get("name"),
            "sort_order": item.get("sortOrder"),
            "input_type_id": item.get("inputType"),
            "mandatory": item.get("mandatory"),
        }
        field = existing.get(item.get("id"))
        if field is not None:
            for key, value in values.items():
                setattr(field, key, value)
        elif not item.get("id"):
            db.session.add(
                models.Field(
                    section_id=section.id,
                    template_instance_id=new_field_instance_id,
                    **values,
                )
            )

    db.session.commit()


def delete_instance_section(
    template_instance: models.WorkupInstance, section: models.Section
) -> None:
    """
    Deletes a section an instance added, with its fields.

    Raises:
        SectionInUseError: If any of its fields has saved values. Saved values are
            kept for the audit trail, so the step should be hidden instead.
        ValueError: If the section doesn't belong to this instance.

    Args:
        template_instance: The instance that owns the section.
        section: The section to delete.
    """
    if section.template_instance_id != template_instance.id:
        raise ValueError("Section does not belong to this instance")
    if _field_has_responses([field.id for field in section.fields]):
        raise SectionInUseError(
            "This step has saved values, so it can't be deleted. Hide it instead."
        )
    db.session.delete(section)
    db.session.commit()
