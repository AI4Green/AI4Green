from sources import models, services


def add(person, workgroup, old_role, new_role):
    """
    Add an entry to the DataAccessChanges table to record a change in data access.
    """

    models.DataAccessChanges.create(
        person_id=person.id,
        workgroup_id=workgroup.id,
        old_role=old_role,
        new_role=new_role,
    )
