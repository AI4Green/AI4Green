import json
from typing import Any, Optional
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
            topic (Optional[str]): The topic to send the message to.
            msg (Optional[Any]): The message to send to the queue.
        """
        self.producer.send(topic, value=msg)
        self.producer.flush()
