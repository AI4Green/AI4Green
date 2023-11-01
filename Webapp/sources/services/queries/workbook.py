from datetime import datetime, timedelta
from typing import List

import pytz
from sources import models
from sources.extensions import db


def list_all() -> List[models.WorkBook]:
    """
    Gets a list of all workbooks

    Returns:
         List of all workbooks
    """
    return (db.session.query(models.WorkBook)).all()
