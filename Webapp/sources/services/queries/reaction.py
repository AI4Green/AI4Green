from datetime import datetime, timedelta
from typing import List

import pytz
from sources import models
from sources.extensions import db


def get_recent_reactions() -> List[models.Reaction]:
    """
    Gets a list of reactions created in the past 28 days

    Returns:
         List of all reactions from past 28 days
    """
    cut_off_date = datetime.now(pytz.timezone("Europe/London")).replace(
        tzinfo=None
    ) - timedelta(days=28)
    return (
        db.session.query(models.Reaction).filter(
            models.Reaction.time_of_creation > cut_off_date
        )
    ).all()


def get_reaction_count() -> int:
    """
    Gets the number of reactions in the database

    Returns:
        The number of reactions in the database
    """
    return db.session.query(models.Reaction).count()
