from sources import models
from sources.extensions import db


def add(person, workgroup, old_role, new_role, workbook=None):
    """
    Add an entry to the DataAccessChanges table to record a change in data access.
    """

    models.DataAccessChanges.create(
        person_id=person.id,
        workgroup_id=workgroup.id,
        workbook_id=workbook.id if workbook else None,
        old_role=old_role,
        new_role=new_role,
    )


def get_changes_from_person(person_id: int):
    """
    Returns the role history for a Person for all WorkGroups
    """
    return (
        db.session.query(models.DataAccessChanges)
        .filter(models.DataAccessChanges.person_id == person_id)
        .all()
    )


def get_changes_from_person_and_workgroup(person_id: int, workgroup_id: int):
    """
    Returns the role history for a Person for a particular WorkGroup
    """
    return (
        db.session.query(models.DataAccessChanges)
        .filter(models.DataAccessChanges.person_id == person_id)
        .filter(models.DataAccessChanges.workgroup_id == workgroup_id)
        .all()
    )


def get_changes_from_workgroup(workgroup_id: int):
    """
    Returns the role history for all Persons for a particular WorkGroup
    """
    return (
        db.session.query(models.DataAccessChanges)
        .filter(models.DataAccessChanges.workgroup_id == workgroup_id)
        .all()
    )
