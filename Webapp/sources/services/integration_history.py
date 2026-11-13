"""
Audit log messages for changes to a workgroup's external integrations, such as
connecting, replacing or removing a ChemInventory API token.

Sent to the "integration-history" queue. The audit_log_saver worker saves the
queues listed in its QUEUE_NAMES setting, so add this queue there.
"""

from dataclasses import asdict, dataclass
from typing import Optional

from flask import current_app, json
from sources.services.message_queue.serde import MessageSerdeMixin

TOPIC = "integration-history"


@dataclass
class IntegrationChangeMessage(MessageSerdeMixin):
    """A change to a workgroup's integration, for the integration-history topic.

    Never include secrets: `token_hint` is only the last few characters of a
    token, as shown to users.
    """

    full_name: str
    email: str
    workgroup: int
    integration: str
    action: str
    date: str
    token_hint: Optional[str] = None

    def serialise(self) -> str:
        """Convert a message into a JSON string.

        Returns:
            str: The message class formatted as a JSON string.
        """
        return json.dumps(asdict(self))

    @staticmethod
    def deserialise(data: dict) -> "IntegrationChangeMessage":
        """Convert a message from a JSON object to a message object.

        Args:
            data (dict): The JSON to decode.

        Returns:
            IntegrationChangeMessage: The message from the JSON.
        """
        return IntegrationChangeMessage(**data)


def send_message(message: IntegrationChangeMessage) -> None:
    """Send a message to the message queue's integration-history topic.

    Args:
        message (IntegrationChangeMessage): The message to send.
    """
    producer = current_app.config["MESSAGE_QUEUE_PRODUCER"]
    producer.send(TOPIC, message.serialise())
