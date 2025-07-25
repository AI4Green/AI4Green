from dataclasses import asdict, dataclass

from flask import current_app, json
from sources.services.message_queue import MessageSerdeMixin


@dataclass
class DataExportMessage(MessageSerdeMixin):
    """Class for creating a kafka message for the data_export_history topic"""

    person: int
    workgroup: int
    workbook: int
    reactions: list[int]
    date: str

    def serialise(self):
        """Convert a message into a JSON object with a schema and payload.

        The schema is an object that lists the fields and their types.
        The payload is an object representing the message class in JSON format.

        Returns:
            dict: The message class formatted with schema and payload.
        """
        schema = {
            "type": "struct",
            "optional": False,
            "fields": [
                {"field": "person", "type": "int32"},
                {"field": "workgroup", "type": "int32"},
                {"field": "workbook", "type": "int32"},
                {"field": "reactions", "type": "list"},
                {"field": "date", "type": "string"},
            ],
        }
        payload = asdict(self)
        serialised = json.dumps({"schema": schema, "payload": payload})
        return serialised

    @staticmethod
    def deserialise(data):
        """Convert a message from a JSON object to a message object.

        The message object is a class decorated with `@dataclass`.

        Args:
            data (dict): The JSON to decode.

        Returns:
            DataExportMessage: The `DataExportMessage` from the JSON.
        """
        return DataExportMessage(**data)


def send_message(message: DataExportMessage):
    """Send a message to the kafka producer in the data_export_history topic.
    Args:
        message (DataExportMessage): The message to send to the queue in the DataExportMessage format
    """
    # producer = current_app.config["MESSAGE_QUEUE_PRODUCER"]
    # producer.send("data_export_history", message.serialise())
