from datetime import datetime
from flask import jsonify, request, send_file
from flask_login import login_required

from sources.services.audit_logs import get_audit_logs, make_log_stream
from sources.decorators import principal_investigator_required

from . import audit_log_bp


@audit_log_bp.route("/download", methods=["GET"])
@login_required
@principal_investigator_required
def download_audit_logs():
    # get arguments form GET parameters
    topic = request.args.get("workbook")
    if topic is None:
        return jsonify({"error": "topic is a required field"}), 400
    workbook = request.args.get("workbook")
    try:
        workbook = int(workbook) if workbook is not None else None
    except (TypeError, ValueError):
        return jsonify({"error": "workbook must an integer, e.g. 1, 2, 3..."}), 400
    start_date = request.args.get("start_date")
    end_date = request.args.get("end_date")

    # retrieve and deserialise the audit logs
    logs = get_audit_logs(
        topic=topic, workbook=workbook, start_date=start_date, end_date=end_date
    )

    # convert the logs into a ZIP file stream
    # file name will be {topic}-{current time}.zip
    # e.g. reaction_editing_history-202506161054.zip
    file_name = f"{topic}-{datetime.now().strftime('%Y%m%d%M%S')}.zip"
    log_stream = make_log_stream(logs=logs, file_name=file_name)

    return send_file(
        log_stream,
        mimetype="application/zip",
        as_attachment=True,
        download_name=file_name,
    )
