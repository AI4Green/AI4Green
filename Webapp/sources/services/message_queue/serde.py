from typing import Any


class MessageSerdeMixin:
    """A mixin for converting classes to JSON serialisable `dict`s that can
    be sent to message queues.
    """

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
