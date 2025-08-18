from dataclasses import asdict, dataclass
import json
from serde import MessageSerdeMixin


@dataclass
class ReactionEditMessage(MessageSerdeMixin):
    """Class for creating a message for the reaction_editing_history topic"""

    full_name: str
    email: str
    workgroup: int
    workbook: int
    reaction: int
    field_name: str
    change_details: dict
    date: str

    def serialise(self) -> str:
        """Convert a message into a JSON string.

        Returns:
            str: The message class formated as a JSON string.
        """
        payload = asdict(self)
        serialised = json.dumps(payload)
        return serialised

    @staticmethod
    def deserialise(data):
        """Convert a message from a JSON object to a message object.

        The message object is a class decorated with `@dataclass`.

        Args:
            data (dict): The JSON to decode.

        Returns:
            ReactionEditMessage: The `ReactionEditMessage` from the JSON.
        """
        return ReactionEditMessage(
            full_name=data["full_name"],
            email=data["email"],
            workgroup=data["workgroup"],
            workbook=data["workbook"],
            reaction=data["reaction"],
            field_name=data["field_name"],
            date=data["date"],
            change_details=data["change_details"],
        )
