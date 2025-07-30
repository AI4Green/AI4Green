import logging
import os
from typing import List
from azure.storage.queue import QueueServiceClient, QueueMessage


def get_messages(sc: QueueServiceClient, queue_name: str) -> List[QueueMessage]:
    """Get messages from Azure Queue Storage with the given queue name.

    Args:
        sc (QueueServiceClient): The service client for Azure Queue Storage.
        queue_name (str): The name of the queue to read from.

    Returns:
        List[QueueMessage]: The list of messages retrieved from the queue.
    """
    return []


def write_messages_to_db(messages: List[QueueMessage], db_connection_string: str):
    """Write the messages from Azure Queue Storage to the database.

    Args:
        messages (List[QueueMessage]): The list of messages to write to the database.
        db_connection_string (str): The connection string for the database.
    """
    pass


def main():
    # default string for testing purposes
    default_str = (
        "DefaultEndpointsProtocol=http;"
        "AccountName=devstoreaccount1;"
        "AccountKey=Eby8vdM02xNOcqFlqUwJPLlmEtlCDXJ1OUzFT50uSRZ6IFsuFq2UVErCz4I6tq/K1SZFPTOtr/KBHBeksoGMGw==;"
        "QueueEndpoint=http://127.0.0.1:10001/devstoreaccount1;"
    )
    connection_str = os.getenv("QUEUE_CONNECTION_STRING", default_str)
    service_client = QueueServiceClient.from_connection_string(connection_str)
    queue_names = os.getenv("QUEUE_NANES", "test")
    queues = queue_names.split(",")

    db_connection_str = os.getenv(
        "DB_CONNECTION_STRING",
        "postgresql://postgres:postgres@localhost:5433/ai4gauditlogtest",
    )

    running = True
    while running:
        try:
            for queue in queues:
                messages = get_messages(service_client, queue)
                write_messages_to_db(messages, db_connection_str)
        except KeyboardInterrupt:
            running = False
            logging.info("Exitiing...")
        except Exception as e:
            running = False
            logging.error(f"An error occurred...{os.linesep}{e}")


if __name__ == "__main__":
    main()
