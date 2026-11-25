"""Sections of a template and their fields."""

from typing import Any, Dict, List, Optional

from sources import models
from sources.extensions import db


class SectionInUseError(Exception):
    """Raised when removing a field that already has saved values."""


def create_section(
    name: str,
    sort_order: Optional[int],
    template: models.Template,
    commit: bool = True,
) -> models.Section:
    """
    Creates a section on a template.

    Args:
        name: Section name.
        sort_order: Position of the section.
        template: The owning template.
        commit: Whether to commit the session.

    Returns:
        The new section.
    """
    section = models.Section(name=name, sort_order=sort_order, template_id=template.id)
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
            updated.append(section)
        else:
            updated.append(
                create_section(
                    name=item.get("name"),
                    sort_order=item.get("sortOrder"),
                    template=template,
                    commit=False,
                )
            )

    db.session.commit()
    return updated


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


def save_fields(section: models.Section, fields_data: List[Dict[str, Any]]) -> None:
    """
    Replaces the fields of a section with the submitted list.

    Fields missing from the list are deleted, fields with an id are updated and
    fields without one are created. Ids that aren't this section's fields are ignored.

    Args:
        section: The section being edited.
        fields_data: Items with "name", "sortOrder", "inputType" (InputType id),
            "mandatory" and, for existing fields, "id".

    Raises:
        SectionInUseError: When a field to delete has saved values in any instance.
    """
    existing = {field.id: field for field in section.fields}
    incoming_ids = {item.get("id") for item in fields_data if item.get("id")}

    ids_to_delete = [field_id for field_id in existing if field_id not in incoming_ids]
    # saved values are kept for the audit trail, so a field with any can't go
    if _field_has_responses(ids_to_delete):
        raise SectionInUseError("This field has saved values, so it can't be deleted.")
    for field_id in ids_to_delete:
        db.session.delete(existing[field_id])

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
            db.session.add(models.Field(section_id=section.id, **values))

    db.session.commit()
