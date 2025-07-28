import json
import time
from collections import defaultdict
from typing import Any

from kafka import KafkaConsumer
from sources.services.message_queue.base import BaseQueueProducer
from sources import services


class MessageSerdeMixin:
    def serialise(self) -> dict:
        """Convert a message into a JSON object with a schema and payload.

        The schema is an object that lists the fields and their types.
        The payload is an object representing the message class in JSON format.

        Raises:
            NotImplementedError: You must implement this method in classes that implement this mixin.

        Returns:
            dict: The message class formated with shcema and payload.
        """
        raise NotImplementedError

    @staticmethod
    def deserialise(data: dict) -> Any:
        """Convert a message from a JSON object to a message object.

        The message object is a class decorated with `@dataclass`.

        Args:
            data (dict): The JSON to decode.

        Raises:
            NotImplementedError: You must implement this method in classes that implement this mixin.

        Returns:
            Any: The message class from the JSON.
        """
        raise NotImplementedError


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

    def poll(self, poll_duration_sec=60):
        """Poll messages from topics and returns them.

        Args:
            poll_duration_sec (int): Length of time to poll kafka for, in seconds
        """
        print(f"Polling messages from topic: {self.topic}")
        messages = []

        end_time = time.time() + poll_duration_sec
        while time.time() < end_time:
            msg_pack = self.consumer.poll(timeout_ms=1000)
            for tp, msgs in msg_pack.items():
                for message in msgs:
                    messages.append(message.value)

        print(f"Collected {len(messages)} messages.")

        return messages


class ReactionEditHistoryProcessor:
    """Processes reaction editing history messages.
    Returns compressed messages to kafka reaction_editing_history_compressed topic.
    """

    def __init__(self, producer: BaseQueueProducer, produce_topic: str):
        self.producer = producer
        self.produce_topic = produce_topic

    def process_and_publish(self, messages):
        if not messages:
            print("No messages to process.")
            return

        print(f"Processing {len(messages)} messages...")

        message_dicts = [json.loads(msg) for msg in messages]

        # group messages by (reaction, field_name, person)
        grouped = defaultdict(list)
        for item in message_dicts:
            key = (item["reaction"], item["field_name"], item["person"])
            grouped[key].append(item)

        # Merge diffs within each group
        for key, messages in grouped.items():
            reaction, field_name, person = key

            # Extract all change_details from messages
            change_details_list = [msg["change_details"] for msg in messages]

            # initially set result to the first change
            change_details_merged = change_details_list[0]

            # loop through subsequent changes (if any) to identify net change
            if len(change_details_list) > 1:
                for diff in change_details_list[1:]:
                    change_details_merged = merge_diffs(change_details_merged, diff)

            if not change_details_merged:  # no net change
                continue

            # put results in original message format
            message = services.reaction_editing_history.ReactionEditMessage(
                person,
                messages[0]["workgroup"],
                messages[0]["workbook"],
                reaction,
                field_name,
                change_details_merged,
                messages[0]["date"],
            )

            # send back to kafka
            self.producer.send(self.produce_topic, message.serialise())


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
