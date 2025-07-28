from typing import Any, Optional


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
