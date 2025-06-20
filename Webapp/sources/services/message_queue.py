import json
import time
from collections import defaultdict
from dataclasses import asdict
from typing import Any, Optional

from flask import current_app
from kafka import KafkaConsumer, KafkaProducer
from sources import services


class BaseQueueProducer:
    def send(self, topic: Optional[str] = None, msg: Optional[Any] = None):
        """Send a message to the message queue with the given topic.

        Args:
            topic (Optional[str]): The topic to send the message to.
            msg (Optional[Any]): The message to send to the queue.
        """
        raise NotImplementedError


class QueueProducer(BaseQueueProducer):
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

    def send(self, topic: Optional[str] = None, msg: Optional[Any] = None):
        """Send a message to the message queue with the given topic.

        Args:
            topic (str, optional): The topic to send the message to.
            msg (Any, optional): The message to send to the queue.
        """
        self.producer.send(topic, value=msg)
        self.producer.flush()


class LoggingQueueProducer(BaseQueueProducer):
    """
    A queue producer designed to be used in testing situations where
    actually sending messages to a message queue is not required.

    This producer simply logs the messages with the Flask logger instance.
    """

    def send(self, topic: Optional[str] = None, msg: Optional[Any] = None):
        """Send a message to the standard output using the built-in logger
        in Flask.

        Args:
            topic (str, optional): This is ignored in this class. Defaults to None.
            msg (Any, optional): The message to send to the queue. Defaults to None.
        """
        current_app.logger.info(msg=msg)


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

    def poll(self, poll_duration_sec=10):  # TODO: adjust poll duration
        """Poll messages from topics and process them.
        Returns compressed messages to kafka reaction_editing_history_compressed topic.

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
    """Processes reaction editing history messages."""

    def __init__(self, producer: BaseQueueProducer):
        self.producer = producer

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
            self.producer.send(
                "reaction_editing_history_compressed", json.dumps(asdict(message))
            )


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
