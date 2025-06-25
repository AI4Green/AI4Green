import os

from dotenv import load_dotenv

# Load environment variables from a .env file if present
load_dotenv()

# Kafka config
KAFKA_HOSTNAME = os.getenv("KAFKA_HOSTNAME", "localhost:9092")
CONSUME_TOPIC = os.getenv("CONSUME_TOPIC", "reaction_editing_history")
PRODUCE_TOPIC = os.getenv("PRODUCE_TOPIC", "reaction_editing_history_compressed")

# Polling config
POLL_INTERVAL_MINS = int(os.getenv("POLL_INTERVAL_MINS", 60))
POLL_DURATION_SEC = int(os.getenv("POLL_DURATION_SEC", 60))
