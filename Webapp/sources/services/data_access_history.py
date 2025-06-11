from dataclasses import asdict, dataclass
from typing import Optional

from flask import current_app, json

from sources.services.message_queue import MessageSerialiserMixin


@dataclass
class DataAccessMessage(MessageSerialiserMixin):
    """Class for creating a kafka message for the data_access_history topic"""

    person: int
    workgroup: int
    old_role: str
    new_role: str
    date: str
    workbook: Optional[int] = None

    def serialise(self):
        """Convert a message into a JSON object with a schema and payload.

        The schema is an object that lists the fields and their types.
        The payload is an object representing the message class in JSON format.

        Returns:
            dict: The message class formated with shcema and payload.
        """
        schema = {
            "type": "struct",
            "optional": False,
            "fields": [
                {"field": "person", "type": "int32"},
                {"field": "workgroup", "type": "int32"},
                {"field": "old_role", "type": "string"},
                {"field": "new_role", "type": "string"},
                {"field": "date", "type": "string"},
                {"field": "workbook", "type": "int32", "optional": True},
            ],
        }
        payload = asdict(self)
        serialised = json.dumps({"schema": schema, "payload": payload})
        return serialised


def send_message(message: DataAccessMessage):
    """Send a message to the kafka producer in the data_access_history topic.
    Args:
        message (DataAccessMessage): The message to send to the queue in the DataAccessMessage format
    """
    producer = current_app.config["MESSAGE_QUEUE_PRODUCER"]
    producer.send("data_access_history", json.dumps(asdict(message)))


# TODO: record account deletion
