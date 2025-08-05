from dataclasses import asdict, dataclass
from typing import Optional

from flask import current_app, json
from sources.services.message_queue.serde import MessageSerdeMixin


@dataclass
class DataAccessMessage(MessageSerdeMixin):
    """Class for creating a kafka message for the data_access_history topic"""

    full_name: str
    email: str
    workgroup: int
    old_role: str
    new_role: str
    date: str
    workbook: Optional[int] = None

    def serialise(self) -> str:
        """Convert a message into a JSON string.

        Returns:
            str: The message class formated as a JSON string.
        """
        payload = asdict(self)
        serialised = json.dumps(payload)
        return serialised

    @staticmethod
    def deserialise(data):
        """Convert a message from a JSON object to a message object.

        The message object is a class decorated with `@dataclass`.

        Args:
            data (dict): The JSON to decode.

        Returns:
            DataAccessMessage: The `DataAccessMessage` from the JSON.
        """
        return DataAccessMessage(**data)


def send_message(message: DataAccessMessage):
    """Send a message to the kafka producer in the data_access_history topic.
    Args:
        message (DataAccessMessage): The message to send to the queue in the DataAccessMessage format
    """
    producer = current_app.config["MESSAGE_QUEUE_PRODUCER"]
    # topic name can only be letters and numbers
    producer.send("dataaccesshistory", message.serialise())
