import json
import time
from collections import defaultdict

from kafka import KafkaConsumer, KafkaProducer


class QueueProducer:
    """This class is a service which is used to send messages to a kafka cluster."""

    def __init__(self, hostname: str):
        """Create a producer which can send messages to a kafka cluster.

        Args:
            hostname (str): The hostname of the cluster, e.g. my.kafka.cluster:9092
        """
        self.producer = KafkaProducer(
            bootstrap_servers=hostname,
            # We'll be sending messages as JSON objects,
            # but kafka accepts messages as binary strings.
            # This lambda converts dicts to strings and then to binary.
            value_serializer=lambda x: json.dumps(x).encode(),
            # client_id="daniel", TODO: find a good value for the client id
        )

    def send(self, topic: str, msg: dict):
        """Send a message to the message queue with the given topic.

        Args:
            topic (str): The topic to send the message to.
            msg (dict): The message to send to the queue.
        """
        self.producer.send(topic, value=msg)
        self.producer.flush()


class QueueConsumer:
    """This class is a service which is used to receive messages from a kafka cluster."""

    def __init__(self, hostname, topic):
        """Create a consumer which can receive messages from a kafka cluster.

        Args:
            hostname (str): The hostname of the cluster, e.g. my.kafka.cluster:9092
            topic (str): The name of the topic
        """
        self.consumer = KafkaConsumer(
            topic,
            bootstrap_servers=hostname,
            group_id="audit-log-consumer-group",
            value_deserializer=lambda x: json.loads(x.decode()),
            auto_offset_reset="earliest",
            enable_auto_commit=True,
        )
        self.topic = topic

    def poll_and_process(self, poll_duration_sec=10):  # TODO: adjust poll duration
        """Poll messages from topics and process them.

        Args:
            poll_duration_sec (int): Length of time to poll kafka for, in seconds
        """
        print(f"Polling messages from topics: {self.topic}")
        messages = []

        end_time = time.time() + poll_duration_sec
        while time.time() < end_time:
            msg_pack = self.consumer.poll(timeout_ms=1000)
            for tp, msgs in msg_pack.items():
                for message in msgs:
                    messages.append(message.value)

        print(f"Collected {len(messages)} messages.")

        # compress messages and send to storage
        processed_messages = process_batch(messages)
        store_batch(processed_messages)


def process_batch(messages):
    """Process a batch of kafka message to find net change in change_details dict.
    Args:
        messages (list): list of json messages from kafka consumer
    """
    print(f"Processing {len(messages)} messages...")

    message_dicts = [json.loads(msg) for msg in messages]
    print(message_dicts)

    # group messages by (reaction, field_name, person)
    grouped = defaultdict(list)
    for item in message_dicts:
        key = (item["reaction"], item["field_name"], item["person"])
        grouped[key].append(item)

    # Merge diffs within each group
    merged_results = []
    for key, messages in grouped.items():
        reaction, field_name, person = key

        # Extract all change_details from messages and load if needed
        change_details_list = []
        for msg in messages:
            change_details = msg["change_details"]
            if isinstance(change_details, str):  # TODO: probably always a string?
                change_details = json.loads(change_details)
            change_details_list.append(change_details)

        # initially set result to the first change
        change_details_merged = change_details_list[0]

        # loop through subsequent changes (if any) to identify net change
        if len(change_details_list) > 1:
            for diff in change_details_list[1:]:
                change_details_merged = merge_diffs(change_details_merged, diff)
        print(change_details_merged)

        # put results in original message format
        results_dict = {
            "person": person,
            "workbook": messages[0]["workbook"],
            "reaction": reaction,
            "field_name": field_name,
            "change_details": change_details_merged,
        }

        merged_results.append(results_dict)

    return merged_results


def store_batch(messages):
    """Send the processed list of messages to storage"""
    return


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
            merged[key] = {
                "old_value": diff1[key]["old_value"],
                "new_value": diff2[key]["new_value"],
            }
        elif key in diff1:
            merged[key] = diff1[key]
        elif key in diff2:
            merged[key] = diff2[key]

    return merged
