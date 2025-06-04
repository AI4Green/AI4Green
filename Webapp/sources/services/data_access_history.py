from dataclasses import asdict, dataclass
from typing import Optional

from flask import current_app, json
from sources import models
from sources.extensions import db


@dataclass
class DataAccessMessage:
    """Class for creating kafka message"""

    person: int
    workgroup: int
    old_role: str
    new_role: str
    workbook: Optional[int] = None


def send_message(message):
    producer = current_app.config["MESSAGE_QUEUE_PRODUCER"]
    producer.send("data_access_history", json.dumps(asdict(message)))


# TODO: record account deletion

'''
def get_history_from_person(person_id: int):
    """
    Returns the role history for a Person for all WorkGroups
    """
    return (
        db.session.query(models.DataAccessHistory)
        .filter(models.DataAccessHistory.person_id == person_id)
        .all()
    )


def get_history_from_workgroup(workgroup_id: int):
    """
    Returns the role history for all Persons for a particular WorkGroup
    """
    return (
        db.session.query(models.DataAccessHistory)
        .filter(models.DataAccessHistory.workgroup_id == workgroup_id)
        .all()
    )
'''
