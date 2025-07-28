from apscheduler.schedulers.blocking import BlockingScheduler
from config import (
    CONSUME_TOPIC,
    KAFKA_HOSTNAME,
    POLL_DURATION_SEC,
    POLL_INTERVAL_MINS,
    PRODUCE_TOPIC,
)
from sources.services.message_queue.consumers import QueueConsumer
from sources.services.message_queue.producers import QueueProducer

from .reaction_edit_history_processor import ReactionEditHistoryProcessor


def collect_and_process_messages(
    consumer_service, poll_duration_sec, processor_service
):
    """Poll Kafka, collect messages, process, and return to Kafka."""
    print("Starting hourly message collection...")
    messages = consumer_service.poll(poll_duration_sec)
    processor_service.process_and_publish(messages)
    print("Finished processing messages.")


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


# set up kafka consumer
consumer = QueueConsumer(
    hostname=KAFKA_HOSTNAME,
    topic=CONSUME_TOPIC,
)

# set up kafka producer and processor
producer = QueueProducer(
    hostname=KAFKA_HOSTNAME,
)
processor = ReactionEditHistoryProcessor(producer, PRODUCE_TOPIC)

# Set up the scheduler
scheduler = BlockingScheduler()
scheduler.add_job(
    collect_and_process_messages,
    args=[consumer, POLL_DURATION_SEC, processor],
    trigger="interval",
    minutes=POLL_INTERVAL_MINS,
)

# Start the scheduler
if __name__ == "__main__":
    print("Starting Kafka consumer scheduler...")
    scheduler.start()
