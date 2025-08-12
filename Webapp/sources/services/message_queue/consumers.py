import json
import time
from typing import List
from azure.storage.queue import QueueServiceClient, QueueMessage
from kafka import KafkaConsumer


class BaseQueueConsumer:

    def poll(self, poll_duration_sec: int = 60):
        """Poll messages from topics and returns them.

        Args:
            poll_duration_sec (int): Length of time to poll for, in seconds
        """
        raise NotImplementedError


class QueueConsumer(BaseQueueConsumer):
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


class AzureQueueConsumer(BaseQueueConsumer):
    """
    A queue consumer for reading messages from Azure Queue Storage.
    """

    def __init__(self, connection_str: str, queue_name: str, max_messages: int = 32):
        """Create a consumer for Azure Queue Storage

        Args:
            connection_str (str): The connection string for the queue storage resource.
            queue_name (str): The name of the queue to read from.
        """
        self.consumer = QueueServiceClient.from_connection_string(connection_str)
        self.queue_name = queue_name
        self.max_messages = max_messages

    def poll(self, poll_duration_sec=60):
        """Poll messages from topics and returns them.

        Args:
            poll_duration_sec (int): Length of time to poll kafka for, in seconds
        """
        messages = []

        client = self.consumer.get_queue_client(self.queue_name)

        end_time = time.time() + poll_duration_sec
        while time.time() < end_time:
            # Read the messages from the queue
            msg_pack = [
                message
                for message in client.receive_messages(max_messages=self.max_messages)
            ]
            for msg in msg_pack:
                messages.append(msg.content)
            # Clear the old mesages
            self._clear_messages(client, msg_pack)

        client.close()

        return messages

    def _clear_messages(self, sc: QueueServiceClient, messages: List[QueueMessage]):
        """Delete messages from the Azure Queue Storage.

        Args:
            sc (QueueServiceClient): The service client for Azure Queue Storage.
            messages (List[QueueMessage]): The list of messages to delete.
        """
        client = sc.get_queue_client(self.queue_name)
        for message in messages:
            client.delete_message(message)
