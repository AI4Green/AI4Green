import logging

from apscheduler.schedulers.blocking import BlockingScheduler
from config import (
    CONSUME_TOPIC,
    POLL_DURATION_SEC,
    POLL_INTERVAL_MINS,
    PRODUCE_TOPIC,
    QUEUE_CONNECTION_STRING,
)
from reaction_edit_history_processor import ReactionEditHistoryProcessor
from message_queue.consumers import AzureQueueConsumer
from message_queue.producers import AzureQueueProducer


def collect_and_process_messages(
    consumer_service, poll_duration_sec, processor_service, logger
):
    """Poll, collect messages, process, and return to message queue."""
    logger.info("Starting hourly message collection...")
    messages = consumer_service.poll(poll_duration_sec)
    processor_service.process_and_publish(messages)
    logger.info("Finished processing messages.")


def main():
    logger = logging.Logger("audit_log_compressor", level=logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # set up consumer
    consumer = AzureQueueConsumer(
        connection_str=QUEUE_CONNECTION_STRING,
        queue_name=CONSUME_TOPIC,
        max_messages=None,
    )

    # set up producer and processor
    producer = AzureQueueProducer(connection_str=QUEUE_CONNECTION_STRING, logger=logger)
    processor = ReactionEditHistoryProcessor(producer, PRODUCE_TOPIC, logger)

    # Set up the scheduler
    scheduler = BlockingScheduler()
    scheduler.add_job(
        collect_and_process_messages,
        args=[consumer, POLL_DURATION_SEC, processor, logger],
        trigger="interval",
        minutes=POLL_INTERVAL_MINS,
    )
    logger.info("Starting consumer scheduler...")
    scheduler.start()


# Set up and start the compressor
if __name__ == "__main__":
    main()
