from typing import Any


class MessageSerdeMixin:
    """A mixin for converting classes to JSON serialisable `dict`s that can
    be sent to message queues.
    """

    def serialise(self) -> str:
        """Convert a message into a JSON string.

        Raises:
            NotImplementedError: You must implement this method in classes that implement this mixin.

        Returns:
            str: The message class formatted as a JSON string.
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
