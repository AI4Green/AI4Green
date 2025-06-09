from Webapp.sources.services.message_queue import QueueConsumer
from flask import current_app
import os
from apscheduler.schedulers.blocking import BlockingScheduler


def collect_and_process_messages():
    """Poll Kafka, collect messages, process, and upload to Minio."""
    print("Starting hourly message collection...")
    consumer.poll_and_process()
    print("Finished processing messages.")

    pass


# set up kafka consumer
consumer = QueueConsumer(
        #hostname=current_app.config["MESSAGE_QUEUE_CONFIG"],
        hostname=os.getenv("MESSAGE_QUEUE_HOSTNAME", "localhost:9092"),
        topic="reaction_editing_history"
)

# Set up the scheduler
scheduler = BlockingScheduler()
scheduler.add_job(collect_and_process_messages, 'interval', seconds=30)

# Start the scheduler
if __name__ == "__main__":
    print("Starting Kafka consumer scheduler...")
    scheduler.start()
