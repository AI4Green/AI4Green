from datetime import datetime, timedelta
from typing import List

import pytz
from sources import models
from sources.extensions import db


def get_recent_workbooks() -> List[models.WorkBook]:
    """
    Gets a list of workbooks created in the past 28 days

    Returns:
         List of all workbooks from past 28 days
    """
    cut_off_date = datetime.now(pytz.timezone("Europe/London")).replace(
        tzinfo=None
    ) - timedelta(days=28)
    return (
        db.session.query(models.WorkBook).filter(
            models.WorkBook.time_of_creation > cut_off_date
        )
    ).all()
