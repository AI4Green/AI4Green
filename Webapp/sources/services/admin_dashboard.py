from datetime import datetime, timedelta
from typing import List

import pytz
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


def get_all_workgroups() -> List[models.WorkGroup]:
    """
    Gets a list of all workgroups in the database

    Returns:
         List of all workgroups
    """
    return db.session.query(models.WorkGroup).all()


def get_all_users() -> List[models.User]:
    """
    Gets a list of all users in the database

    Returns:
         List of all users
    """
    return db.session.query(models.User).all()


def get_new_workgroup_requests() -> List[models.WorkGroupRequest]:
    """
    Gets a list of all workgroup requests in the database

    Returns:
         List of all workgroup requests
    """
    return (
        db.session.query(models.WorkGroupRequest)
        .filter(models.WorkGroupRequest.status == "active")
        .all()
    )


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


def get_reaction_count() -> int:
    """
    Gets the number of reactions in the database

    Returns:
        The number of reactions in the database
    """
    return db.session.query(models.Reaction).count()


def get_compound_count() -> int:
    """
    Gets the number of compounds in the database

    Returns:
        The number of compounds in the database
    """
    return db.session.query(models.Compound).count()
