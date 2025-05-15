from sources import models
from sources.extensions import db


def add(person, workgroup, old_role, new_role, workbook=None):
    """
    Add an entry to the DataAccessHistory table to record a change in data access.
    """

    models.DataAccessHistory.create(
        person_id=person.id,
        workgroup_id=workgroup.id,
        workbook_id=workbook.id if workbook else None,
        old_role=old_role,
        new_role=new_role,
    )


def get_history_from_person(person_id: int):
    """
    Returns the role history for a Person for all WorkGroups
    """
    return (
        db.session.query(models.DataAccessHistory)
        .filter(models.DataAccessHistory.person_id == person_id)
        .all()
    )


def get_history_from_person_and_workgroup(person_id: int, workgroup_id: int):
    """
    Returns the role history for a Person for a particular WorkGroup
    """
    return (
        db.session.query(models.DataAccessHistory)
        .filter(models.DataAccessHistory.person_id == person_id)
        .filter(models.DataAccessHistory.workgroup_id == workgroup_id)
        .all()
    )


def get_history_from_workgroup(workgroup_id: int):
    """
    Returns the role history for all Persons for a particular WorkGroup
    """
    return (
        db.session.query(models.DataAccessHistory)
        .filter(models.DataAccessHistory.workgroup_id == workgroup_id)
        .all()
    )
