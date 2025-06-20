import os

from apscheduler.schedulers.blocking import BlockingScheduler

from Webapp.sources.services.message_queue import (
    QueueConsumer,
    QueueProducer,
    ReactionEditHistoryProcessor,
)


def collect_and_process_messages():
    """Poll Kafka, collect messages, process, and return to Kafka."""
    print("Starting hourly message collection...")
    messages = consumer.poll()
    processor.process_and_publish(messages)
    print("Finished processing messages.")

    pass


# set up kafka consumer
consumer = QueueConsumer(
    hostname=os.getenv("MESSAGE_QUEUE_HOSTNAME", "localhost:9092"),
    topic="reaction_editing_history",
)

# set up kafka producer and processor
producer = QueueProducer(
    hostname=os.getenv("MESSAGE_QUEUE_HOSTNAME", "localhost:9092"),
)
processor = ReactionEditHistoryProcessor(producer)

# Set up the scheduler
scheduler = BlockingScheduler()
scheduler.add_job(collect_and_process_messages, "interval", seconds=30)

# Start the scheduler
if __name__ == "__main__":
    print("Starting Kafka consumer scheduler...")
    scheduler.start()
