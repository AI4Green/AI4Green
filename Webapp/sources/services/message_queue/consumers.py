import json
import time
from kafka import KafkaConsumer


class QueueConsumer:
    """This class is a service which is used to receive messages from a kafka cluster."""

    def __init__(self, hostname, topic):
        """Create a consumer which can receive messages from a kafka cluster.

        Args:
            hostname (str): The hostname of the cluster, e.g. my.kafka.cluster:9092
            topic (str): The name of the topic
        """
        self.consumer = KafkaConsumer(
            topic,
            bootstrap_servers=hostname,
            group_id="audit-log-consumer-group",
            value_deserializer=lambda x: json.loads(x.decode()),
            auto_offset_reset="earliest",
            enable_auto_commit=True,
        )
        self.topic = topic

    def poll(self, poll_duration_sec=60):
        """Poll messages from topics and returns them.

        Args:
            poll_duration_sec (int): Length of time to poll kafka for, in seconds
        """
        print(f"Polling messages from topic: {self.topic}")
        messages = []

        end_time = time.time() + poll_duration_sec
        while time.time() < end_time:
            msg_pack = self.consumer.poll(timeout_ms=1000)
            for tp, msgs in msg_pack.items():
                for message in msgs:
                    messages.append(message.value)

        print(f"Collected {len(messages)} messages.")

        return messages
