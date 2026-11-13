"""
Read-only integration with ChemInventory (https://www.cheminventory.net).

A workgroup can store a ChemInventory API token. The inventory page can then
search the workgroup's ChemInventory by CAS number or name, and show where the
matching containers are, grouped by location. Only searches are made: loading
the whole inventory needs a call per location, which is far too slow for
large inventories.

Tokens:
- Stored encrypted with Fernet, bound to their workgroup, using the keys in
  CHEMINVENTORY_ENCRYPTION_KEY (in Azure, an App Setting referencing Key
  Vault). The setting may list several comma-separated keys for rotation: the
  first encrypts, and tokens under older keys are re-encrypted when next used.
- Never logged, and never sent back to the browser. Only the last few
  characters are shown, so people can tell which token is in use.
- Tokens for any ChemInventory account type are accepted, but a token for a
  read-only user is recommended: an administrator's token can change or
  delete anything in the inventory if it were ever leaked. Saving an
  administrator's token is logged as a warning.
- Each ChemInventory token acts as one ChemInventory user, so every member of
  the workgroup sees what that user can see.

A search is two calls (the inventory's id, then the search itself), with an
overall time limit.

ChemInventory API: every call is a POST of JSON with `authtoken` in the body,
and every response is `{"status": "success" | "error", "data": ...}`.
"""

import logging
import re
import time
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

import requests
from cryptography.fernet import Fernet, InvalidToken, MultiFernet
from flask import current_app
from flask_login import current_user
from sources import models, services
from sources.extensions import db
from sqlalchemy.exc import IntegrityError

logger = logging.getLogger(__name__)

# Time limit for a single call, and for a whole search (up to two calls)
REQUEST_TIMEOUT_SECONDS = 20
SEARCH_TIME_LIMIT_SECONDS = 30

SEARCH_TYPES = ("cas", "name")
# CAS registry numbers: 2-7 digits, 2 digits, 1 check digit (as the client's
# CAS_PATTERN in components/inventory/utils.jsx)
CAS_PATTERN = re.compile(r"[1-9]\d{1,6}-\d{2}-\d", re.ASCII)
MIN_NAME_LENGTH = 2
MAX_QUERY_LENGTH = 100
# Containers kept from one search; a common solvent can match hundreds
MAX_SEARCH_RESULTS = 200

TOKEN_HINT_LENGTH = 4
MAX_TOKEN_LENGTH = 512
MAX_LOGGED_TEXT_LENGTH = 200

# ChemInventory account types (from /general/getdetails) with more permissions
# than this integration needs: 2 = Group Administrator, 3 = Order Handler.
# (4 = User and 6 = Read-only are enough; the reaction integration will need a
# standard user, to record usage.)
PRIVILEGED_ACCOUNT_TYPES = (2, 3)

NOT_AVAILABLE_MESSAGE = "ChemInventory isn't set up on this server."
NOT_CONNECTED_MESSAGE = "This workgroup hasn't connected a ChemInventory account."
UNREACHABLE_MESSAGE = "ChemInventory couldn't be reached. Please try again later."
TOO_SLOW_MESSAGE = "ChemInventory took too long to respond. Please try again later."
UNEXPECTED_RESPONSE_MESSAGE = "ChemInventory sent a response we couldn't read."
REJECTED_MESSAGE = "ChemInventory didn't accept the API token."
UNREADABLE_TOKEN_MESSAGE = (
    "The stored ChemInventory token can't be read. Please connect ChemInventory again."
)


class ChemInventoryError(Exception):
    """A problem talking to ChemInventory. The message is safe to show users."""


class ChemInventoryRejectedError(ChemInventoryError):
    """ChemInventory answered, but refused the request (e.g. a bad token)."""


class ChemInventoryNotConnectedError(Exception):
    """The workgroup has no ChemInventory connection."""


class ChemInventoryUnavailableError(Exception):
    """The integration is switched off or misconfigured on this server."""


class ChemInventoryValidationError(ValueError):
    """Invalid input from the user. The message is safe to show users."""


def _loggable(value: Any) -> str:
    """
    Text from outside the app, made safe to log: shortened, and quoted so
    line breaks can't forge extra log lines.

    Args:
        value: The text.

    Returns:
        The shortened, quoted text.
    """
    return repr(str(value)[:MAX_LOGGED_TEXT_LENGTH])


# --- Token encryption --------------------------------------------------------


def _encryption_keys() -> List[Fernet]:
    """
    The configured encryption keys, newest first.

    Returns:
        A cipher per key.

    Raises:
        ChemInventoryUnavailableError: If there are no keys, or one isn't a
            valid Fernet key.
    """
    setting = current_app.config.get("CHEMINVENTORY_ENCRYPTION_KEY") or ""
    keys = [key.strip() for key in setting.split(",") if key.strip()]
    if not keys:
        raise ChemInventoryUnavailableError(NOT_AVAILABLE_MESSAGE)
    try:
        return [Fernet(key) for key in keys]
    except ValueError:
        logger.error("CHEMINVENTORY_ENCRYPTION_KEY contains an invalid Fernet key")
        raise ChemInventoryUnavailableError(NOT_AVAILABLE_MESSAGE)


def _api_url() -> str:
    """
    The ChemInventory API's base URL.

    Returns:
        The URL, without a trailing slash.

    Raises:
        ChemInventoryUnavailableError: If it isn't HTTPS, as tokens must never
            be sent unencrypted.
    """
    url = current_app.config.get("CHEMINVENTORY_API_URL") or ""
    if not url.lower().startswith("https://"):
        logger.error("CHEMINVENTORY_API_URL must be an https:// URL")
        raise ChemInventoryUnavailableError(NOT_AVAILABLE_MESSAGE)
    return url.rstrip("/")


def is_available() -> bool:
    """
    Whether the integration is switched on and correctly configured.

    Returns:
        True if there are valid encryption keys and an HTTPS API URL.
    """
    try:
        _encryption_keys()
        _api_url()
    except ChemInventoryUnavailableError:
        return False
    return True


def encrypt_token(workgroup_id: int, token: str) -> str:
    """
    Encrypts an API token for storage, with the newest key. The workgroup id is
    encrypted with it, so a stored token can't be moved to another workgroup.

    Args:
        workgroup_id: The workgroup the token belongs to.
        token: The plain token.

    Returns:
        The encrypted token, as text.
    """
    newest_key = _encryption_keys()[0]
    return newest_key.encrypt(f"{workgroup_id}:{token}".encode()).decode()  # noqa: E231


def decrypt_token(workgroup_id: int, encrypted_token: str) -> str:
    """
    Decrypts a stored API token, with whichever configured key encrypted it.

    Args:
        workgroup_id: The workgroup the token is stored for.
        encrypted_token: The token as stored.

    Returns:
        The plain token.

    Raises:
        ChemInventoryError: If it can't be decrypted (e.g. its key has been
            removed), or it was encrypted for a different workgroup.
    """
    try:
        plaintext = MultiFernet(_encryption_keys()).decrypt(encrypted_token.encode())
    except InvalidToken:
        raise ChemInventoryError(UNREADABLE_TOKEN_MESSAGE)

    stored_workgroup_id, separator, token = plaintext.decode().partition(":")
    if not separator or stored_workgroup_id != str(workgroup_id):
        logger.error("A stored ChemInventory token belongs to a different workgroup")
        raise ChemInventoryError(UNREADABLE_TOKEN_MESSAGE)
    return token


def _uses_old_key(encrypted_token: str) -> bool:
    """
    Whether a stored token was encrypted with a key other than the newest.

    Args:
        encrypted_token: The token as stored.

    Returns:
        True if it should be re-encrypted with the newest key.
    """
    try:
        _encryption_keys()[0].decrypt(encrypted_token.encode())
    except InvalidToken:
        return True
    return False


# --- HTTP client -------------------------------------------------------------


@dataclass(frozen=True)
class _Client:
    """
    Calls the ChemInventory API with one token, until a deadline.
    """

    base_url: str
    deadline: float  # time.monotonic() value after which no call is started
    repr = False  # keeps the token out of logs and tracebacks that print the client
    token: str = field(repr=False)

    def post(self, path: str, body: Optional[Dict[str, Any]] = None) -> Any:
        """
        Calls an API endpoint.

        Args:
            path: The endpoint path, e.g. "/location/load".
            body: Request parameters, other than the token.

        Returns:
            The response's `data`.

        Raises:
            ChemInventoryRejectedError: If ChemInventory returns an error.
            ChemInventoryError: If it can't be reached in time, or the response
                can't be read.
        """
        remaining = self.deadline - time.monotonic()
        if remaining <= 0:
            raise ChemInventoryError(TOO_SLOW_MESSAGE)

        try:
            response = requests.post(
                f"{self.base_url}{path}",
                json={**(body or {}), "authtoken": self.token},
                timeout=min(REQUEST_TIMEOUT_SECONDS, remaining),
                # a redirect could resend the body, with the token, elsewhere
                allow_redirects=False,
            )
        except requests.Timeout:
            logger.warning("ChemInventory request to %s timed out", path)
            raise ChemInventoryError(TOO_SLOW_MESSAGE)
        except requests.RequestException as error:
            # the exception text includes the URL but not the request body
            logger.warning(
                "ChemInventory request to %s failed: %s", path, _loggable(error)
            )
            raise ChemInventoryError(UNREACHABLE_MESSAGE)

        if 300 <= response.status_code < 400:
            logger.warning(
                "ChemInventory %s redirected (HTTP %s), which isn't followed",
                path,
                response.status_code,
            )
            raise ChemInventoryError(UNEXPECTED_RESPONSE_MESSAGE)

        try:
            payload = response.json()
        except ValueError:
            logger.warning(
                "ChemInventory %s returned non-JSON (HTTP %s)",
                path,
                response.status_code,
            )
            raise ChemInventoryError(UNEXPECTED_RESPONSE_MESSAGE)

        status = payload.get("status") if isinstance(payload, dict) else None
        if status == "error":
            # ChemInventory's error text is logged rather than shown, as it may
            # describe the account
            logger.info(
                "ChemInventory %s returned an error: %s",
                path,
                _loggable(payload.get("data")),
            )
            raise ChemInventoryRejectedError(REJECTED_MESSAGE)
        if status != "success" or not response.ok:
            logger.warning(
                "ChemInventory %s returned an unexpected response (HTTP %s)",
                path,
                response.status_code,
            )
            raise ChemInventoryError(UNEXPECTED_RESPONSE_MESSAGE)
        return payload.get("data")


def _client(token: str, time_limit: float) -> _Client:
    """
    Builds a client for a token, using the configured API URL.

    Args:
        token: The plain API token.
        time_limit: Seconds from now after which no further call is started.

    Returns:
        The client.
    """
    return _Client(
        base_url=_api_url(), deadline=time.monotonic() + time_limit, token=token
    )


# --- Connection --------------------------------------------------------------


def _read_token(value: Any) -> str:
    """
    Validates an API token from the request body.

    Args:
        value: The submitted token.

    Returns:
        The token, trimmed.

    Raises:
        ChemInventoryValidationError: If it's missing or malformed.
    """
    if not isinstance(value, str) or not value.strip():
        raise ChemInventoryValidationError("Enter an API token.")
    token = value.strip()
    if len(token) > MAX_TOKEN_LENGTH or any(c.isspace() for c in token):
        raise ChemInventoryValidationError("That doesn't look like an API token.")
    return token


def _as_int(value: Any) -> Optional[int]:
    """A ChemInventory number (which may be a string) as an int, or None."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _truncate(value: Any, max_length: int) -> Optional[str]:
    """A value as trimmed text no longer than max_length, or None if empty."""
    if value is None:
        return None
    text = str(value).strip()
    return text[:max_length] or None


def _is_set(value: Any) -> bool:
    """Whether a ChemInventory flag (1, "1", true or "true") is set."""
    return str(value).strip().lower() in ("1", "true")


def _is_privileged(user: Dict[str, Any]) -> bool:
    """
    Whether a ChemInventory user has more permissions than this integration
    needs: an administrator or order handler.

    Args:
        user: The `user` from /general/getdetails.

    Returns:
        True for administrators and order handlers.
    """
    return (
        _as_int(user.get("accounttype")) in PRIVILEGED_ACCOUNT_TYPES
        or _is_set(user.get("siteadmin"))
        or _is_set(user.get("enterpriseadmin"))
    )


def check_token(token: str) -> Dict[str, Optional[str]]:
    """
    Checks a token works, by asking ChemInventory who it belongs to.

    Args:
        token: The plain API token.

    Returns:
        `user` (the ChemInventory user's name) and `inventory` (the inventory
        name), either of which may be None.

    Raises:
        ChemInventoryRejectedError: If ChemInventory refuses the token.
        ChemInventoryError: If ChemInventory can't be reached.
    """
    details = _client(token, REQUEST_TIMEOUT_SECONDS).post("/general/getdetails")
    details = details if isinstance(details, dict) else {}
    user = details.get("user") if isinstance(details.get("user"), dict) else {}
    inventory = (
        details.get("inventory") if isinstance(details.get("inventory"), dict) else {}
    )

    if _is_privileged(user):
        logger.warning(
            "A ChemInventory token for an administrator or order handler is "
            "being saved; a read-only user's token is recommended"
        )

    full_name = " ".join(
        part for part in (user.get("firstname"), user.get("lastname")) if part
    )
    return {
        "user": _truncate(full_name, 255),
        "inventory": _truncate(inventory.get("name") or user.get("inventoryname"), 255),
    }


def get_connection(
    workgroup: models.WorkGroup,
) -> Optional[models.ChemInventoryConnection]:
    """
    Gets the workgroup's ChemInventory connection.

    Args:
        workgroup: The workgroup.

    Returns:
        The connection, or None if it hasn't connected.
    """
    return (
        db.session.query(models.ChemInventoryConnection)
        .filter(models.ChemInventoryConnection.workgroup_id == workgroup.id)
        .one_or_none()
    )


def get_connection_status(workgroup: models.WorkGroup) -> Dict[str, Any]:
    """
    The workgroup's connection as the settings modal shows it. Never includes
    the token.

    Args:
        workgroup: The workgroup.

    Returns:
        `available` (switched on for this server), `configured`, `canManage`
        (the user may change the connection), and for a configured
        workgroup `tokenHint`, `cheminventoryUser`, `inventoryName`,
        `connectedBy` and `connectedAt`.
    """
    connection = get_connection(workgroup)
    status: Dict[str, Any] = {
        "available": is_available(),
        "configured": connection is not None,
        "canManage": services.inventory.can_edit(workgroup),
    }
    if connection is not None:
        status.update(
            {
                "tokenHint": connection.token_hint,
                "cheminventoryUser": connection.cheminventory_user,
                "inventoryName": connection.inventory_name,
                "connectedBy": connection.created_by.fullname
                if connection.created_by
                else None,
                "connectedAt": connection.updated_at.isoformat()
                if connection.updated_at
                else None,
            }
        )
    return status


def _record_change(
    workgroup: models.WorkGroup, action: str, token_hint: Optional[str]
) -> None:
    """
    Records a change to the workgroup's ChemInventory token in the audit log.

    Args:
        workgroup: The workgroup.
        action: "connected", "token_replaced" or "disconnected".
        token_hint: The last few characters of the token concerned.
    """
    message = services.integration_history.IntegrationChangeMessage(
        full_name=current_user.fullname,
        email=current_user.email,
        workgroup=workgroup.id,
        integration="cheminventory",
        action=action,
        token_hint=token_hint,
        date=datetime.now(timezone.utc).isoformat(),
    )
    services.integration_history.send_message(message)


def _apply_token(
    connection: models.ChemInventoryConnection,
    encrypted_token: str,
    token_hint: str,
    details: Dict[str, Optional[str]],
) -> None:
    """Sets a connection's token and the details ChemInventory gave for it."""
    connection.encrypted_token = encrypted_token
    connection.token_hint = token_hint
    connection.cheminventory_user = details["user"]
    connection.inventory_name = details["inventory"]
    connection.created_by_id = current_user.id
    connection.updated_at = datetime.now(timezone.utc)


def save_connection(
    workgroup: models.WorkGroup, data: Dict[str, Any]
) -> models.ChemInventoryConnection:
    """
    Connects the workgroup to ChemInventory, or replaces its token. The token
    is checked with ChemInventory before it's stored.

    Args:
        workgroup: The workgroup.
        data: The request body, with `token`.

    Returns:
        The saved connection.

    Raises:
        ChemInventoryValidationError: If the token is missing or malformed.
        ChemInventoryRejectedError: If ChemInventory refuses the token.
        ChemInventoryError: If ChemInventory can't be reached.
        ChemInventoryUnavailableError: If the integration is switched off.
    """
    token = _read_token(data.get("token"))
    encrypted_token = encrypt_token(workgroup.id, token)  # fails early if off
    details = check_token(token)
    token_hint = token[-TOKEN_HINT_LENGTH:]

    connection = get_connection(workgroup)
    action = "token_replaced" if connection else "connected"
    if connection is None:
        connection = models.ChemInventoryConnection(workgroup_id=workgroup.id)
        db.session.add(connection)
    _apply_token(connection, encrypted_token, token_hint, details)
    try:
        db.session.commit()
    except IntegrityError:
        # someone else connected this workgroup at the same moment: replace
        # their token with this one
        db.session.rollback()
        connection = get_connection(workgroup)
        if connection is None:
            raise
        action = "token_replaced"
        _apply_token(connection, encrypted_token, token_hint, details)
        db.session.commit()

    _record_change(workgroup, action, token_hint)
    return connection


def delete_connection(workgroup: models.WorkGroup) -> None:
    """
    Disconnects the workgroup from ChemInventory, deleting its stored token.

    Args:
        workgroup: The workgroup.

    Raises:
        ChemInventoryNotConnectedError: If it isn't connected.
    """
    connection = get_connection(workgroup)
    if connection is None:
        raise ChemInventoryNotConnectedError(NOT_CONNECTED_MESSAGE)
    token_hint = connection.token_hint
    db.session.delete(connection)
    db.session.commit()
    _record_change(workgroup, "disconnected", token_hint)


# --- Search ------------------------------------------------------------------


def _inventory_id(client: _Client) -> int:
    """
    The id of the ChemInventory inventory a token searches. Searches need it,
    and only /general/getdetails gives it.

    Args:
        client: A client for the workgroup's token.

    Returns:
        The inventory id.

    Raises:
        ChemInventoryError: If ChemInventory can't be reached, or doesn't say
            which inventory the token belongs to.
    """
    details = client.post("/general/getdetails")
    details = details if isinstance(details, dict) else {}
    inventory = details.get("inventory")
    user = details.get("user")
    inventory_id = _as_int(inventory.get("id")) if isinstance(inventory, dict) else None
    if inventory_id is None and isinstance(user, dict):
        # the user's active inventory, in case the inventory object is missing
        inventory_id = _as_int(user.get("inventory"))
    if inventory_id is None:
        logger.warning("ChemInventory /general/getdetails gave no inventory id")
        raise ChemInventoryError(UNEXPECTED_RESPONSE_MESSAGE)
    return inventory_id


def _number(value: Any) -> Optional[float]:
    """A ChemInventory number (which may be a string) as a float, or None."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _iso_date(value: Any) -> Optional[str]:
    """
    A ChemInventory date as YYYY-MM-DD, or None if it isn't in that format.
    Other formats (ChemInventory can use the user's chosen date format) are
    dropped rather than guessed at.
    """
    text = _truncate(value, 10)
    if not text:
        return None
    try:
        return date.fromisoformat(text).isoformat()
    except ValueError:
        return None


def _read_search(data: Dict[str, Any]) -> Tuple[str, str]:
    """
    Validates a search from the request body.

    Args:
        data: The request body, with `type` ("cas" or "name") and `query`.

    Returns:
        The search type and the trimmed query.

    Raises:
        ChemInventoryValidationError: If either is missing or malformed.
    """
    search_type = data.get("type")
    if search_type not in SEARCH_TYPES:
        raise ChemInventoryValidationError("Choose to search by CAS number or name.")

    query = data.get("query")
    query = query.strip() if isinstance(query, str) else ""
    if search_type == "cas":
        if not CAS_PATTERN.fullmatch(query):
            raise ChemInventoryValidationError("Enter a CAS number like 67-64-1.")
    elif not MIN_NAME_LENGTH <= len(query) <= MAX_QUERY_LENGTH:
        raise ChemInventoryValidationError(
            f"Enter a name of {MIN_NAME_LENGTH} to {MAX_QUERY_LENGTH} characters."
        )
    return search_type, query


def _serialise_container(raw: Dict[str, Any]) -> Dict[str, Any]:
    """
    A container from a ChemInventory search, with only the fields the page
    shows. All of ChemInventory's container fields are read here, so if its
    format changes (it can depend on the inventory's table layout), this is the
    one place to update. Free-text comments are left out, as they may hold
    notes not meant for every member.

    Args:
        raw: A container from /search/execute.

    Returns:
        The container.
    """
    # `name` is the container's name (usually the substance's); `substance` is
    # only a numeric substance ID, so it isn't shown
    return {
        "id": raw.get("id"),
        "name": _truncate(raw.get("name"), 255) or "Unnamed substance",
        "cas": _truncate(raw.get("cas"), 50),
        "quantity": _number(raw.get("size")),
        "unit": _truncate(raw.get("unit"), 16) or "",
        "barcode": _truncate(raw.get("barcode"), 100),
        "supplier": _truncate(raw.get("supplier"), 255),
        "dateReceived": _iso_date(raw.get("dateacquired")),
    }


def _location_path(location_text: str) -> List[str]:
    """
    Splits ChemInventory's location text into its parts.

    Args:
        location_text: The full path, e.g. "Cupboards > Cupboard 1".

    Returns:
        The non-empty parts, outermost first, e.g. ["Cupboards", "Cupboard 1"],
        or the whole text as the only part if it has none.
    """
    parts = [part.strip() for part in location_text.split(">")]
    return [part for part in parts if part] or [location_text]


def _group_by_location(raw_containers: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Groups search results into one entry per location, for the location tiles.

    Args:
        raw_containers: Containers from /search/execute.

    Returns:
        Locations sorted by name, each with `id`, `name` (the full path, e.g.
        "Cupboards > Cupboard 1"), `path` (the same path as a list, e.g.
        ["Cupboards", "Cupboard 1"]) and its `containers`.
    """
    locations: Dict[str, Dict[str, Any]] = {}
    for raw in raw_containers:
        location_id = _as_int(raw.get("location"))
        name = _truncate(raw.get("locationtext"), 255) or "Unknown location"
        key = (
            f"id:{location_id}"  # noqa: E231
            if location_id is not None
            else f"name:{name}"  # noqa: E231
        )
        if key not in locations:
            locations[key] = {
                "id": location_id,
                "name": name,
                "path": _location_path(name),
                "containers": [],
            }
        locations[key]["containers"].append(_serialise_container(raw))
    return sorted(locations.values(), key=lambda loc: loc["name"].lower())


def search(workgroup: models.WorkGroup, data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Searches the workgroup's ChemInventory for containers of a substance, by
    CAS number or name (ChemInventory includes synonyms in name searches).

    Args:
        workgroup: The workgroup.
        data: The request body, with `type` ("cas" or "name") and `query`.

    Returns:
        `type`, `query`, `locations` (one per location holding a match, each
        with its `containers`), `containerCount`, and `truncated`, which is
        True when there were more than MAX_SEARCH_RESULTS matches and only
        the first were kept.

    Raises:
        ChemInventoryValidationError: If the search is missing or malformed.
        ChemInventoryUnavailableError: If the integration is switched off.
        ChemInventoryNotConnectedError: If the workgroup hasn't connected.
        ChemInventoryError: If ChemInventory can't be reached, or refuses.
    """
    search_type, query = _read_search(data)
    if not is_available():
        raise ChemInventoryUnavailableError(NOT_AVAILABLE_MESSAGE)
    connection = get_connection(workgroup)
    if connection is None:
        raise ChemInventoryNotConnectedError(NOT_CONNECTED_MESSAGE)

    token = decrypt_token(workgroup.id, connection.encrypted_token)
    if _uses_old_key(connection.encrypted_token):
        connection.encrypted_token = encrypt_token(workgroup.id, token)
        db.session.commit()

    client = _client(token, SEARCH_TIME_LIMIT_SECONDS)
    result = client.post(
        "/search/execute",
        {
            "inventory": _inventory_id(client),
            "type": search_type,
            "contents": query,
        },
    )
    containers = result.get("containers") if isinstance(result, dict) else None
    containers = [c for c in containers or [] if isinstance(c, dict)]
    kept = containers[:MAX_SEARCH_RESULTS]

    return {
        "type": search_type,
        "query": query,
        "locations": _group_by_location(kept),
        "containerCount": len(kept),
        "truncated": len(containers) > MAX_SEARCH_RESULTS,
    }
