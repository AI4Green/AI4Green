import time
from typing import List
from azure.storage.queue import QueueClient, QueueServiceClient, QueueMessage


class BaseQueueConsumer:

    def poll(self, poll_duration_sec: int = 60):
        """Poll messages from topics and returns them.

        Args:
            poll_duration_sec (int): Length of time to poll for, in seconds
        """
        raise NotImplementedError


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

    def _clear_messages(self, client: QueueClient, messages: List[QueueMessage]):
        """Delete messages from the Azure Queue Storage.

        Args:
            client (QueueClient): The client for Azure Queue Storage.
            messages (List[QueueMessage]): The list of messages to delete.
        """
        for message in messages:
            client.delete_message(message)
