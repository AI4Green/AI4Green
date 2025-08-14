import os

from dotenv import load_dotenv

# Load environment variables from a .env file if present
load_dotenv()

# message queue config
QUEUE_CONNECTION_STRING = os.getenv("QUEUE_CONNECTION_STRING", "localhost:9092")
CONSUME_TOPIC = os.getenv("CONSUME_TOPIC", "reaction-editing-history")
PRODUCE_TOPIC = os.getenv("PRODUCE_TOPIC", "reaction-editing-history-compressed")

# Polling config
POLL_INTERVAL_MINS = int(os.getenv("POLL_INTERVAL_MINS", 60))
POLL_DURATION_SEC = int(os.getenv("POLL_DURATION_SEC", 60))
