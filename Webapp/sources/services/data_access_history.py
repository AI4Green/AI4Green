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
    workbook: Optional[int] = None


def send_message(message):
    # producer = current_app.config["MESSAGE_QUEUE_PRODUCER"]
    # producer.send("data_access_history", json.dumps(asdict(message)))
    return


# TODO: record account deletion
