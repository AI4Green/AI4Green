from dataclasses import asdict, dataclass
from typing import Optional

from flask import current_app, json


@dataclass
class DataAccessMessage:
    """Class for creating a kafka message for the data_access_history topic"""

    person: int
    workgroup: int
    old_role: str
    new_role: str
    date: str
    workbook: Optional[int] = None


def send_message(message: DataAccessMessage):
    """Send a message to the kafka producer in the data_access_history topic.
    Args:
        message (DataAccessMessage): The message to send to the queue in the DataAccessMessage format
    """
    producer = current_app.config["MESSAGE_QUEUE_PRODUCER"]
    producer.send("data_access_history", json.dumps(asdict(message)))


# TODO: record account deletion
