import os

from apscheduler.schedulers.blocking import BlockingScheduler
from config import (
    CONSUME_TOPIC,
    KAFKA_HOSTNAME,
    POLL_DURATION_SEC,
    POLL_INTERVAL_MINS,
    PRODUCE_TOPIC,
)

from Webapp.sources.services.message_queue import (
    QueueConsumer,
    QueueProducer,
    ReactionEditHistoryProcessor,
)


def collect_and_process_messages(
    consumer_service, poll_duration_sec, processor_service
):
    """Poll Kafka, collect messages, process, and return to Kafka."""
    print("Starting hourly message collection...")
    messages = consumer_service.poll(poll_duration_sec)
    processor_service.process_and_publish(messages)
    print("Finished processing messages.")


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
