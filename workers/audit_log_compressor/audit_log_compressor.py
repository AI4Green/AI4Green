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
from sources.services.message_queue.consumers import AzureQueueConsumer
from sources.services.message_queue.producers import AzureQueueProducer


def collect_and_process_messages(
    consumer_service, poll_duration_sec, processor_service, logger
):
    """Poll, collect messages, process, and return to message queue."""
    logger.info("Starting hourly message collection...")
    messages = consumer_service.poll(poll_duration_sec)
    processor_service.process_and_publish(messages)
    logger.info("Finished processing messages.")


def merge_diffs(diff1, diff2):
    """Returns net change between two change_details dicts
    Args:
        diff1 (dict): must have nested "old_value" and "new_value" keys
        diff2 (dict): must have nested "old_value" and "new_value" keys
    """
    merged = {}

    all_keys = set(diff1.keys()) | set(diff2.keys())

    for key in all_keys:
        if key in diff1 and key in diff2:
            if diff1[key]["old_value"] != diff2[key]["new_value"]:
                merged[key] = {
                    "old_value": diff1[key]["old_value"],
                    "new_value": diff2[key]["new_value"],
                }
        elif key in diff1:
            merged[key] = diff1[key]
        elif key in diff2:
            merged[key] = diff2[key]

    return merged


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
    producer = AzureQueueProducer(
        connection_str=QUEUE_CONNECTION_STRING,
    )
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
