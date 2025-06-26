from datetime import datetime

from flask import current_app, jsonify, request, send_file
from flask_login import login_required
from sources.decorators import principal_investigator_required
from sources.services.audit_logs import (
    get_audit_logs,
    get_human_readable_ids,
    make_log_stream,
)

from . import audit_log_bp


@audit_log_bp.route("/audit_log/<workgroup>/download", methods=["GET"])
@login_required
@principal_investigator_required
def download_audit_logs(workgroup: str):
    if not current_app.config["USE_KAFKA"]:
        return (
            jsonify(
                {
                    "success": False,
                    "error": "Kafka is disabled. Refer to the README to enable it.",
                }
            ),
            400,
        )

    # get arguments form GET parameters
    topic = request.args.get("topic")
    if topic is None:
        return jsonify({"error": "topic is a required field"}), 400
    start_date = request.args.get("start_date")
    end_date = request.args.get("end_date")

    # retrieve and deserialise the audit logs
    logs = get_audit_logs(
        topic=topic, workgroup_name=workgroup, start_date=start_date, end_date=end_date
    )

    # make the names of users, workgroups and workbooks human readable
    get_human_readable_ids(logs=logs)

    # convert the logs into a ZIP file stream
    # file name will be {topic}-{current time}
    # e.g. reaction_editing_history-202506161054
    file_name = f"{topic}-{datetime.now().strftime('%Y%m%d%M%S')}"
    log_stream = make_log_stream(logs=logs, file_name=file_name)

    return send_file(
        log_stream,
        mimetype="application/zip",
        as_attachment=True,
        download_name=file_name,
    )
