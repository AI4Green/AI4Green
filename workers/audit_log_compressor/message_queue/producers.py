from logging import Logger
from typing import Any, Optional

from azure.storage.queue import QueueServiceClient


class BaseQueueProducer:
    """
    The interface for message queue producer classes. It provides a
    `send` method in which the inheriting class defines logic to send
    messages to some kind of message queue service (e.g. RabbitMQ, Redis, etc.).
    """

    def send(self, topic: Optional[str] = None, msg: Optional[Any] = None):
        """Send a message to the message queue with the given topic.

        Args:
            topic (Optional[str]): The topic to send the message to.
            msg (Optional[Any]): The message to send to the queue.
        """
        raise NotImplementedError


class AzureQueueProducer(BaseQueueProducer):
    """
    A queue producer for sending messages to Azure Queue Storage.
    """

    def __init__(self, connection_str: str, logger: Logger):
        """Create a producer for Azure Queue Storage

        Args:
            connection_str (str): The connection string for the queue storage resource.
        """
        self.producer = QueueServiceClient.from_connection_string(connection_str)
        self.logger = logger

    def send(self, topic: Optional[str] = None, msg: Optional[Any] = None):
        """Send a message to Azure Queue Storage.

        **NB**: the arguments topic and msg are *required*. The `Optional` type
        annotation is for compatibility with the base class only.

        Args:
            topic (Optional[str], optional):
                The name of the quueue to send the message. Defaults to None.
            msg (Optional[Any], optional): The message to put in the queue. Defaults to None.

        Raises:
            ValueError: 'topic' is required. This should be the name of the queue.
            ValueError: 'msg' is required.
        """
        if topic is None:
            raise ValueError(
                "'topic' is required. This should be the name of the queue."
            )
        if msg is None:
            raise ValueError("'msg' is required.")

        client = self.producer.get_queue_client(topic)

        # check if topic exists
        queues = self.producer.list_queues(name_starts_with=topic)
        queue_names = [queue["name"] for queue in queues]
        if topic not in queue_names:
            # create the queue if it doesn't already exist
            client.create_queue()

        try:
            client.send_message(msg)
        except Exception as e:
            self.logger.error(str(e))
