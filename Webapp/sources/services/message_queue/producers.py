from typing import Any, Optional

from azure.storage.queue import QueueServiceClient
from kafka import KafkaProducer


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
            value_serializer=lambda x: x.encode(),
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
        from flask import current_app

        current_app.logger.info(msg=msg)


class AzureQueueProducer(BaseQueueProducer):
    """
    A queue producer for sending messages to Azure Queue Storage.
    """

    def __init__(self, connection_str: str):
        """Create a producer for Azure Queue Storage

        Args:
            connection_str (str): The connection string for the queue storage resource.
        """
        self.producer = QueueServiceClient.from_connection_string(connection_str)

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
        for queue in queues:
            if queue["name"] != topic:
                # create the queue if it doesn't already exist
                client.create_queue()

        try:
            client.send_message(msg)
        except Exception as e:
            # import to avoid context based conflicts
            from flask import current_app

            current_app.logger.error(msg=str(e))
