# import io
# import json
# import re
# import zipfile
# from datetime import datetime
# from typing import Any, List, Optional

# from flask import current_app
# from minio import Minio
# from sources import services

# client: Minio = current_app.config["MINIO_CLIENT"]


# def _extract_logs(logs: list) -> List[Any]:
#     """Decode the logs from a S3 storage object.

#     Args:
#         logs (list): The raw list of logs from the S3 API containing the log records.

#     Returns:
#         List[Any]: The log records.
#     """
#     extracted = []
#     for line in logs:
#         extracted.append(json.loads(line))
#     return extracted


# def _filter_names_by_date(
#     item_names: List[str],
#     start_date: Optional[str] = None,
#     end_date: Optional[str] = None,
# ) -> List[str]:
#     """Filter names of log files based on the date in their file paths.

#     If `start_date` and `end_date` are given, only logs between and including thos
#     dates will be selected.

#     If only `start_date` is given, only logs after and including that date will be
#     selected.

#     If only `end_date` is given, only logs before and including that date will be
#     selected.

#     If neither are given, all logs will be returned.

#     Args:
#         item_names (List[str]): The list of log file names to filter.
#         start_date (Optional[str], optional):
#             A date in YYYY-MM-DD format to start filtering from.
#             Defaults to None.
#         end_date (Optional[str], optional):
#             A date in YYYY-MM-DD format to end filtering at.
#             Defaults to None.

#     Returns:
#         List[str]: The names of the log files that passed the filter.
#     """
#     if not start_date and not end_date:
#         # No date filtering needed, return all items
#         return item_names

#     date_format = "%Y-%m-%d"
#     date_pattern = re.compile(r"date=(\d{4}-\d{2}-\d{2})")
#     filtered_items = []

#     for item in item_names:
#         date_match = date_pattern.search(item)
#         if not date_match:
#             continue  # No date=YYYY-MM-DD pattern found

#         try:
#             date_str = date_match.group(1)
#             date_obj = datetime.strptime(date_str, date_format)
#         except (ValueError, IndexError):
#             continue  # skip malformed dates

#         if start_date:
#             start = datetime.strptime(start_date, date_format)
#             if date_obj < start:
#                 continue

#         if end_date:
#             end = datetime.strptime(end_date, date_format)
#             if date_obj > end:
#                 continue

#         filtered_items.append(item)

#     return filtered_items


# def get_audit_logs(
#     topic: str,
#     workgroup_name: Optional[str] = None,
#     start_date: Optional[str] = None,
#     end_date: Optional[str] = None,
# ) -> List[Any]:
#     """Get the audit logs for a specific topic. The logs can be filtered by workgroup and date.

#     If no workgroup is given, logs for all workgroups will be selected.

#     If `start_date` and `end_date` are given, only logs between and including thos
#     dates will be selected.

#     If only `start_date` is given, only logs after and including that date will be
#     selected.

#     If only `end_date` is given, only logs before and including that date will be
#     selected.

#     If neither are given, all logs will be returned.

#     Args:
#         topic (str): The topic to retrieve logs for.
#         workgroup_name (Optional[str], optional):
#             The workgroup name to get the logs for. Defaults to None.
#         start_date (Optional[str], optional):
#             Get logs after and including this date. Defaults to None.
#         end_date (Optional[str], optional):
#             Get logs before and including this date. Defaults to None.

#     Returns:
#         List[Any]: The logs deserialised into their appropriate class.
#     """
#     bucket_name = current_app.config["MINIO_AUDIT_LOG_BUCKET"]

#     workgroup_id = services.workgroup.from_name(workgroup_name).id

#     # Read all the object info in the bucket under the provided topic
#     # If a workgroup is given, filter by workgroup too
#     prefix = (
#         f"topics/{topic}"
#         if workgroup_id is None
#         else f"topics/{topic}/workgroup={workgroup_id}"
#     )
#     log_files = list(client.list_objects(bucket_name, recursive=True, prefix=prefix))
#     log_file_names = [obj.object_name for obj in log_files]

#     # Filter the file names by whether they fit in the date bounds
#     log_file_names = _filter_names_by_date(log_file_names, start_date, end_date)

#     # Get logs from files
#     logs = []
#     for name in log_file_names:
#         response = client.get_object(bucket_name=bucket_name, object_name=name)
#         logs.extend(_extract_logs(response.readlines()))
#         response.close()
#         response.release_conn()

#     return logs


# def make_log_stream(logs: List[Any], file_name: str):
#     """
#     Creates an in-memory ZIP file containing the given JSON data.

#     Args:
#         logs: The JSON-serializable data to include in the ZIP.
#         file_name: The name of the JSON file inside the ZIP archive.

#     Returns:
#         A BytesIO object containing the ZIP file data.
#     """
#     # Prepare in-memory bytes buffer
#     zip_buffer = io.BytesIO()

#     # Create a ZIP file in memory
#     with zipfile.ZipFile(zip_buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
#         # Convert JSON to bytes and write into the ZIP archive
#         json_bytes = json.dumps(logs, indent=2).encode("utf-8")
#         zf.writestr(file_name, json_bytes)

#     # Reset the stream's position to the beginning
#     zip_buffer.seek(0)
#     return zip_buffer


# def get_human_readable_ids(logs: List[dict]):
#     for log in logs:
#         # Look up user, workgroup and workbook
#         person = services.person.from_id(log["person"]).user
#         workgroup = services.workgroup.from_id(log["workgroup"])
#         workbook = services.workbook.get(log["workbook"])

#         # get the human readable names
#         fullname = person.fullname if person is not None else "Deleted User"
#         workgroup_name = workgroup.name if workgroup else "Deleted Workgroup"
#         workbook_name = workbook.name if workbook is not None else None

#         # Alter the records to show the human readable names
#         log["person"] = fullname
#         log["workgroup"] = workgroup_name
#         log["workbook"] = workbook_name
