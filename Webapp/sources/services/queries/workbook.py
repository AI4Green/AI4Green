from datetime import datetime, timedelta
from typing import List

import pytz
from sources import models
from sources.extensions import db


def get_all_workbooks() -> List[models.WorkBook]:
    """
    Gets a list of workbooks created in the past 28 days

    Returns:
         List of all workbooks from past 28 days
    """
    return (db.session.query(models.WorkBook)).all()
