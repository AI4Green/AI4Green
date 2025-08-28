from datetime import datetime

from flask import jsonify, request, send_file
from flask_login import login_required
from sources.decorators import principal_investigator_required
from sources.services.audit_logs import (
    get_audit_logs,
    get_human_readable_ids,
    make_log_stream,
)
from sources.services import workgroup as wg

from . import audit_log_bp


@audit_log_bp.route("/audit_log/<workgroup>/download", methods=["GET"])
@login_required
@principal_investigator_required
def download_audit_logs(workgroup: str):
    # get arguments form GET parameters
    topic = request.args.get("topic")
    if topic is None:
        return jsonify({"error": "topic is a required field"}), 400

    # records are saved with the workgroup id but this function receives
    # the name, so get the ID
    workgroup_id = None
    workgroup_obj = wg.from_name(workgroup)
    if workgroup_obj:
        workgroup_id = workgroup_obj.id

    start_date = request.args.get("start_date")
    # convert to datetime
    start_date = datetime.strptime(start_date, "%Y-%m-%d") if start_date else None

    end_date = request.args.get("end_date")
    # convert to datetime
    end_date = datetime.strptime(end_date, "%Y-%m-%d") if end_date else None
    if end_date is not None:
        end_date = end_date.replace(hour=23, minute=59, second=59)

    # retrieve and deserialise the audit logs
    logs = get_audit_logs(
        topic=topic,
        workgroup_name=str(workgroup_id),  # throws an error if not converted to str
        start_date=start_date,
        end_date=end_date,
    )

    # make the names of users, workgroups and workbooks human readable
    get_human_readable_ids(logs=logs)

    # convert the logs into a ZIP file stream
    # file name will be {topic}-{current time}.json
    # e.g. reaction_editing_history-202506161054.json
    file_name = f"{topic}-{datetime.now().strftime('%Y%m%d%M%S')}.json"
    log_stream = make_log_stream(logs=logs, file_name=file_name)

    return send_file(
        log_stream,
        mimetype="application/zip",
        as_attachment=True,
        # `send_file` creates a file around the stream,
        # so change the extension to make it download properly
        download_name=file_name.replace(".json", ".zip"),
    )
