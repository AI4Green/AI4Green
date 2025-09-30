import datetime as dt
import json
import logging
import os
import time
from typing import List

from azure.storage.queue import QueueMessage, QueueServiceClient
from models import AuditLogEvent
from sqlalchemy import Engine, create_engine
from sqlalchemy.orm import Session


def get_messages(
    sc: QueueServiceClient, queue_name: str, max_messages: int
) -> List[QueueMessage]:
    """Get messages from Azure Queue Storage with the given queue name.

    Args:
        sc (QueueServiceClient): The service client for Azure Queue Storage.
        queue_name (str): The name of the queue to read from.
        max_messages (int, optional): The maximum number of messages to be read from the queue.

    Returns:
        List[QueueMessage]: The list of messages retrieved from the queue.
    """
    client = sc.get_queue_client(queue_name)
    messages = [
        message for message in client.receive_messages(max_messages=max_messages)
    ]
    client.close()
    return messages


def clear_messages(
    sc: QueueServiceClient, queue_name: str, messages: List[QueueMessage]
):
    """Delete messages from the Azure Queue Storage.

    Args:
        sc (QueueServiceClient): The service client for Azure Queue Storage.
        queue_name (str): The name of the queue to delete from.
        messages (List[QueueMessage]): The list of messages to delete.
    """
    client = sc.get_queue_client(queue_name)
    for message in messages:
        client.delete_message(message)
    client.close()


def write_messages_to_db(engine: Engine, messages: List[QueueMessage], event_type: str):
    """Write the messages from Azure Queue Storage to the database.

    Args:
        engine (Engine): The engine that connects to the database.
        messages (List[QueueMessage]): The list of messages to write to the database.
        event_type (str): The type of the event.
    """
    # Create the session
    with Session(engine) as session:
        session.begin()
        try:
            for message in messages:
                event_type = event_type
                # get the content from the queue message
                content = json.loads(message.content)
                # generate the table entry
                log_event = AuditLogEvent(
                    event_type=event_type,
                    event_time=dt.datetime.now(),
                    message=content,
                )
                session.add(log_event)
        except:
            # rollback if anything goes wrong
            session.rollback()
            raise
        else:
            # if all goes well, commit the changes
            session.commit()


def main():
    logger = logging.Logger("audit_log_saver", level=logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # default string for testing purposes
    default_str = (
        "DefaultEndpointsProtocol=http;"
        "AccountName=devstoreaccount1;"
        "AccountKey=Eby8vdM02xNOcqFlqUwJPLlmEtlCDXJ1OUzFT50uSRZ6IFsuFq2UVErCz4I6tq/K1SZFPTOtr/KBHBeksoGMGw==;"
        "QueueEndpoint=http://127.0.0.1:10001/devstoreaccount1;"
    )
    connection_str = os.getenv("QUEUE_CONNECTION_STRING", default_str)
    queue_names = os.getenv("QUEUE_NAMES", "test-queue")
    queues = queue_names.split(",")
    db_connection_str = os.getenv(
        "DB_CONNECTION_STRING",
        "postgresql://postgres:postgres@localhost:5434/ai4gauditlog",
    )
    # Establish DB connection
    engine = create_engine(db_connection_str)

    max_messages = int(os.getenv("MAX_MESSAGES", 32))

    service_client = QueueServiceClient.from_connection_string(connection_str)

    running = True
    logger.info(f"Saving messages from queues: '{queues}'.")
    while running:
        retrieved_messages = 0
        try:
            for queue in queues:
                messages = get_messages(service_client, queue, max_messages)
                retrieved_messages += len(messages)
                write_messages_to_db(engine, messages, queue)
                clear_messages(service_client, queue, messages)
        except KeyboardInterrupt:
            running = False
            logger.info("Exitiing...")
        except Exception as e:
            running = False
            logger.error(f"An error occurred...{os.linesep}{e}")
        # If no messages were retrived from the queues, wait 5 seconds
        # to reduce the number of calls to read the queue per minute
        # to save money
        if retrieved_messages == 0:
            time.sleep(5)

    service_client.close()


if __name__ == "__main__":
    main()
