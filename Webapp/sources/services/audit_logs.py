from datetime import datetime
import json
import re
from typing import Any, Callable, List, Optional
from flask import current_app
from minio import Minio


def _extract_logs(logs: list, deserialiser: Callable) -> List[Any]:
    """Decode the logs from a S3 storage object.

    Args:
        logs (list): The raw list of logs from the S3 API containing the log records.
        deserialiser (Callable):
            The function to deserialise the log records. This should be the `deserialise`
            static method on a message class with the `MessageSerdeMixin`.

    Returns:
        List[Any]: The log records.
    """
    extracted = []
    for line in logs:
        extracted.append(deserialiser(json.loads(line)))
    return extracted


def _filter_names_by_date(
    item_names: List[str],
    start_date: Optional[str] = None,
    end_date: Optional[str] = None,
) -> List[str]:
    """Filter names of log files based on the date in their file paths.

    If `start_date` and `end_date` are given, only logs between and including thos
    dates will be selected.

    If only `start_date` is given, only logs after and including that date will be
    selected.

    If only `end_date` is given, only logs before and including that date will be
    selected.

    If neither are given, all logs will be returned.

    Args:
        item_names (List[str]): The list of log file names to filter.
        start_date (Optional[str], optional):
            A date in YYYY-MM-DD format to start filtering from.
            Defaults to None.
        end_date (Optional[str], optional):
            A date in YYYY-MM-DD format to end filtering at.
            Defaults to None.

    Returns:
        List[str]: The names of the log files that passed the filter.
    """
    date_format = "%Y-%m-%d"
    date_pattern = re.compile(r"date=(\d{4}-\d{2}-\d{2})")
    filtered_items = []

    for item in item_names:
        date_match = date_pattern.search(item)
        if not date_match:
            continue  # No date=YYYY-MM-DD pattern found

        try:
            date_str = date_str = date_match.group(1)
            date_obj = datetime.strptime(date_str, date_format)
        except (ValueError, IndexError):
            continue  # skip malformed dates

        if start_date:
            start = datetime.strptime(start_date, date_format)
            if date_obj < start:
                continue

        if end_date:
            end = datetime.strptime(end_date, date_format)
            if date_obj > end:
                continue

        filtered_items.append(item)

    return filtered_items


def get_audit_logs(
    topic: str,
    deserialiser: Callable,
    workbook: Optional[int] = None,
    start_date: Optional[str] = None,
    end_date: Optional[str] = None,
) -> List[Any]:
    """Get the audit logs for a specific topic. The logs can be filtered by workbook and date.

    If no workbook is given, logs for all workbooks will be selected.

    If `start_date` and `end_date` are given, only logs between and including thos
    dates will be selected.

    If only `start_date` is given, only logs after and including that date will be
    selected.

    If only `end_date` is given, only logs before and including that date will be
    selected.

    If neither are given, all logs will be returned.

    Args:
        topic (str): The topic to retrieve logs for.
        deserialiser (Callable):
            The function to deserialise the log records. This should be the `deserialise`
            static method on a message class with the `MessageSerdeMixin`.
        workbook (Optional[int], optional):
            The workbook to get the logs for. Defaults to None.
        start_date (Optional[str], optional):
            Get logs after and including this date. Defaults to None.
        end_date (Optional[str], optional):
            Get logs before and including this date. Defaults to None.

    Returns:
        List[Any]: The logs deserialised into their appropriate class.
    """
    # Set up the Minio client
    client = Minio(
        current_app.config["MINIO_HOST"],
        access_key=current_app.config["MINIO_ACCESS_KEY"],
        secret_key=current_app.config["MINIO_SECRET_KEY"],
        secure=current_app.config["MINIO_SECURE"],
    )
    bucket_name = current_app.config["MINIO_AUDIT_LOG_BUCKET"]

    # Read all the object info in the bucket under the provided topic
    # If a workbook is given, filter by workbook too
    prefix = (
        f"topics/{topic}" if workbook is None else f"topics/{topic}/workbook={workbook}"
    )
    log_files = list(client.list_objects(bucket_name, recursive=True, prefix=prefix))
    log_file_names = [obj.object_name for obj in log_files]

    # Filter the file names by whether they fit in the date bounds
    log_file_names = _filter_names_by_date(log_file_names, start_date, end_date)

    # Get logs from files
    logs = []
    for name in log_file_names:
        response = client.get_object(bucket_name=bucket_name, object_name=name)
        logs.extend(_extract_logs(response.readlines(), deserialiser=deserialiser))
        response.close()
        response.release_conn()

    return logs
