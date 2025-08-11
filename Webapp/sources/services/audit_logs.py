import io
import json
import zipfile
from datetime import datetime
from sqlalchemy import cast, String, Integer
from typing import Any, List, Optional

from sources import services
from sources.extensions import db
from sources.models.audit_log import AuditLogEvent


def get_audit_logs(
    topic: str,
    workgroup_name: Optional[str] = None,
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
) -> List[Any]:
    """Get the audit logs for a specific topic. The logs can be filtered by workgroup and date.

    If no workgroup is given, logs for all workgroups will be selected.

    If `start_date` and `end_date` are given, only logs between and including thos
    dates will be selected.

    If only `start_date` is given, only logs after and including that date will be
    selected.

    If only `end_date` is given, only logs before and including that date will be
    selected.

    If neither are given, all logs will be returned.

    Args:
        topic (str): The topic to retrieve logs for.
        workgroup_name (Optional[str], optional):
            The workgroup name to get the logs for. Defaults to None.
        start_date (Optional[datetime], optional):
            Get logs after and including this date. Defaults to None.
        end_date (Optional[datetime], optional):
            Get logs before and including this date. Defaults to None.

    Returns:
        List[Any]: The logs deserialised into their appropriate class.
    """

    query = db.session.query(AuditLogEvent).filter(AuditLogEvent.event_type == topic)

    # Filter by workgroup if provided
    if workgroup_name is not None:
        query = query.filter(
            cast(AuditLogEvent.message["workgroup"], String) == workgroup_name
        )

    # Apply date filters
    if start_date is not None:
        query = query.filter(AuditLogEvent.event_time >= start_date)
    if end_date is not None:
        query = query.filter(AuditLogEvent.event_time <= end_date)

    query = query.order_by(AuditLogEvent.event_time.asc())

    # Collect the audit log me
    logs = [log.message for log in query.all()]

    return logs


def make_log_stream(logs: List[Any], file_name: str):
    """
    Creates an in-memory ZIP file containing the given JSON data.

    Args:
        logs: The JSON-serializable data to include in the ZIP.
        file_name: The name of the JSON file inside the ZIP archive.

    Returns:
        A BytesIO object containing the ZIP file data.
    """
    # Prepare in-memory bytes buffer
    zip_buffer = io.BytesIO()

    # Create a ZIP file in memory
    with zipfile.ZipFile(zip_buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        # Convert JSON to bytes and write into the ZIP archive
        json_bytes = json.dumps(logs, indent=2).encode("utf-8")
        zf.writestr(file_name, json_bytes)

    # Reset the stream's position to the beginning
    zip_buffer.seek(0)
    return zip_buffer


def get_human_readable_ids(logs: List[dict]):
    for log in logs:
        # Look up user, workgroup and workbook
        workgroup = services.workgroup.from_id(log["workgroup"])
        workbook = services.workbook.get(log["workbook"])

        # get the human readable names
        workgroup_name = workgroup.name if workgroup else "Deleted Workgroup"
        workbook_name = workbook.name if workbook is not None else None

        # Alter the records to show the human readable names
        log["workgroup"] = workgroup_name
        log["workbook"] = workbook_name
