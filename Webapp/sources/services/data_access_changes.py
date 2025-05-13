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
