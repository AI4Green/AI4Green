import json
from typing import Any, Optional

from flask import current_app
from kafka import KafkaProducer


class BaseQueueProducer:
    def send(self, topic: Optional[str] = None, msg: Optional[Any] = None):
        """Send a message to the message queue with the given topic.

        Args:
            topic (Optional[str]): The topic to send the message to.
            msg (Optional[Any]): The message to send to the queue.
        """
        raise NotImplementedError


class QueueProducer(BaseQueueProducer):
    """This class is a service which is used to send messages to a kafka cluster."""

    def __init__(self, hostname: str):
        """Create a producer which can send messages to a kafka cluster.

        Args:
            hostname (str): The hostname of the cluster, e.g. my.kafka.cluster:9092
        """
        self.producer = KafkaProducer(
            bootstrap_servers=hostname,
            # We'll be sending messages as JSON objects,
            # but kafka accepts messages as binary strings.
            # This lambda converts dicts to strings and then to binary.
            value_serializer=lambda x: json.dumps(x).encode(),
            # client_id="daniel", TODO: find a good value for the client id
        )

    def send(self, topic: Optional[str] = None, msg: Optional[Any] = None):
        """Send a message to the message queue with the given topic.

        Args:
            topic (str, optional): The topic to send the message to.
            msg (Any, optional): The message to send to the queue.
        """
        self.producer.send(topic, value=msg)
        self.producer.flush()


class LoggingQueueProducer(BaseQueueProducer):
    """
    A queue producer designed to be used in testing situations where
    actually sending messages to a message queue is not required.

    This producer simply logs the messages with the Flask logger instance.
    """

    def send(self, topic: Optional[str] = None, msg: Optional[Any] = None):
        """Send a message to the standard output using the built-in logger
        in Flask.

        Args:
            topic (str, optional): This is ignored in this class. Defaults to None.
            msg (Any, optional): The message to send to the queue. Defaults to None.
        """
        current_app.logger.info(msg=msg)


class MessageSerdeMixin:
    def serialise(self) -> dict:
        """Convert a message into a JSON object with a schema and payload.

        The schema is an object that lists the fields and their types.
        The payload is an object representing the message class in JSON format.

        Raises:
            NotImplementedError: You must implement this method in classes that implement this mixin.

        Returns:
            dict: The message class formated with shcema and payload.
        """
        raise NotImplementedError

    @staticmethod
    def deserialise(data: dict) -> Any:
        """Convert a message from a JSON object to a message object.

        The message object is a class decorated with `@dataclass`.

        Args:
            data (dict): The JSON to decode.

        Raises:
            NotImplementedError: You must implement this method in classes that implement this mixin.

        Returns:
            Any: The message class from the JSON.
        """
        raise NotImplementedError
