from typing import List

from sources import models
from sources.extensions import db


def get_compound_data_error_reports() -> List[models.CompoundDataErrorReport]:
    """
    Gets a list of compound database error reports in the database

    Returns:
         List of all compound database error reports
    """
    # Add time cut off if list grows too large
    return db.session.query(models.CompoundDataErrorReport).all()
