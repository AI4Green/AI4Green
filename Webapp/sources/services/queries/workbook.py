from typing import List

from sources import models
from sources.extensions import db


def list_all() -> List[models.WorkBook]:
    """
    Gets a list of all workbooks

    Returns:
         List of all workbooks
    """
    return (db.session.query(models.WorkBook)).all()


def get(primary_key: int) -> models.WorkBook:
    """
    Gets the workbook object from the primary key

    Args:
        primary_key: the integer id

    Returns:
        the workbook with the specified primary key
    """
    return (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.id == primary_key)
        .first()
    )


def get_workbook_from_group_book_name_combination(
    workgroup_name: str, workbook_name: str
) -> models.WorkBook:
    """
    Gets the workbook object from the workbook name and the name of the workgroup it belongs to.
    Args:
        workgroup_name: the name of the workgroup that the workbook belongs to
        workbook_name: the name of the workbook that we are looking for

    Returns:
        The workbook which matches the name and group combination
    """

    return (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.name == workbook_name)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )
