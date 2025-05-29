import json
from kafka import KafkaProducer


class QueueProducer:
    def __init__(self, hostname: str):
        self.producer = KafkaProducer(
            bootstrap_servers=hostname,
            value_serializer=lambda x: json.dumps(x).encode(),
            client_id="daniel",
        )

    def send(self, topic: str, msg: dict):
        """Send a message to the message queue with the given topic.

        Args:
            topic (str): The topic to send the message to.
            msg (dict): The message to send to the queue.
        """
        self.producer.send(topic, value=msg)
        self.producer.flush()
