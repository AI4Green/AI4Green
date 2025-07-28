import json
from collections import defaultdict

from sources.services.message_queue.base import BaseQueueProducer
from sources import services


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
