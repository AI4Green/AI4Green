from typing import List

from sources import models
from sources.extensions import db


def get_workbook_result_types(
    workgroup_name: str, workbook_name: str
) -> List[models.ReactionResultType]:
    """
    Gets reaction result types for a given workgroup and workbook, including default result types.

    Args:
        workgroup_name: The name of the workgroup containing the workbook to query.
        workbook_name: The name of the workbook containing the reaction result types.

    Returns:
        List of ReactionResultType instances, including both workbook-specific values
        and pre-seeded default values (where workgroup/workbook are null).
    """
    return (
        db.session.query(models.ReactionResultType)
        .outerjoin(
            models.WorkGroup,
            models.ReactionResultType.workgroup_id == models.WorkGroup.id,
        )
        .outerjoin(
            models.WorkBook, models.ReactionResultType.workbook_id == models.WorkBook.id
        )
        .filter(
            (models.WorkGroup.name == workgroup_name)
            | models.ReactionResultType.workgroup_id.is_(None),
            (models.WorkBook.name == workbook_name)
            | models.ReactionResultType.workbook_id.is_(None),
        )
        .all()
    )
