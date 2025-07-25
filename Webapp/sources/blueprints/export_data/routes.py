import threading
from datetime import datetime

import sources.services.data_export.requests
from flask import (
    Response,
    current_app,
    flash,
    jsonify,
    redirect,
    render_template,
    request,
    url_for,
)
from flask_login import login_required
from sources import models, services
from sources.auxiliary import (
    get_notification_number,
    get_workgroups,
    security_member_workgroup_workbook,
    security_pi_workgroup,
)

from . import export_data_bp


@export_data_bp.route("/export_data/home", methods=["GET", "POST"])
@login_required
@export_data_bp.doc(security="sessionAuth")
def export_data_home():
    """
    Renders the data export home page

    Returns:
        flask.Response: The rendered template for the data export home page
    """
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "data_export/data_export_home.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@export_data_bp.route("/export_data/new_request", methods=["POST", "GET"])
@login_required
@export_data_bp.doc(security="sessionAuth")
def new_data_export_request():
    """
    Initiates a data export request if user has permission to export from the selected workbook.

    Returns:
        flask.Response: A JSON response with the status of the request
    """
    export_request = services.data_export.requests.NewRequest()
    permission_status = export_request.check_permissions()
    if permission_status == "permission accepted":
        export_request.submit_request()
    return jsonify(permission_status)


@export_data_bp.route("/export_data/request_response/<token>", methods=["GET", "POST"])
@login_required
@export_data_bp.doc(security="sessionAuth")
def export_data_request_response(token: str) -> Response:
    """
    Verify the token and the requestor and either redirect or render the page with the request information

    Args:
        token: The token to verify the request

    Returns:
        flask.Response: The rendered template for the data export request page
    """
    # verify and give user the request data or if they fail verification send them home.
    verification = services.data_export.requests.RequestLinkVerification(token)
    if verification.verify_approver_link() == "failed":
        return redirect(url_for("main.index"))
    # get request workbook names from list to render template message
    workbooks_names = [wb.name for wb in verification.data_export_request.workbooks]
    workbooks = ", ".join(workbooks_names)
    return render_template(
        "data_export/data_export_request.html",
        data_export_request=verification.data_export_request,
        workbooks=workbooks,
        workgroups=get_workgroups(),
        notification_number=get_notification_number(),
    )


@export_data_bp.route("/export_data/export_denied", methods=["POST", "GET"])
@login_required
@export_data_bp.doc(security="sessionAuth")
def export_denied():
    """
    Update the database with the denial of the request.

    Returns:
        flask.Response: A redirect to the home page with a success message
    """
    data_export_request = models.DataExportRequest.query.get(request.json["exportID"])
    request_status = services.data_export.requests.RequestStatus(data_export_request)
    request_status.deny()
    flash("Data export denied!")
    return redirect(url_for("main.index"))


@export_data_bp.route("/export_data/export_approved", methods=["POST", "GET"])
@login_required
@export_data_bp.doc(security="sessionAuth")
def export_approved():
    """
    Update the database with the approval. Then check if all approvers have done so and start export.

    Returns:
        flask.Response: A redirect to the home page with a success message
    """
    data_export_request = models.DataExportRequest.query.get(request.json["exportID"])
    request_status = services.data_export.requests.RequestStatus(data_export_request)
    request_status.accept()
    request_status.update_status()
    if data_export_request.status.value == "APPROVED":
        # initiate data export release with threads then notify requestor after successful data export generation.
        export_process = services.data_export.export.initiate
        task_thread = threading.Thread(
            target=export_process,
            args=[current_app.app_context(), data_export_request.id],
        )

        task_thread.start()

    flash("Data export approved!")
    return redirect(url_for("main.index"))


@export_data_bp.route("/export_data/request_download/<token>", methods=["POST", "GET"])
@login_required
@export_data_bp.doc(security="sessionAuth")
def request_download(token: str) -> Response:
    """
    Verify the token and the requestor and either redirect or render the page with the link to download the export file.

    Args:
        token: The token to verify the request

    Returns:
        flask.Response: The rendered template for the data export download page with a link to download the export.
    """
    # verify and give user the request data or if they fail verification send them home.
    verification = services.data_export.requests.RequestLinkVerification(token)
    if verification.verify_requestor_link() == "failed":
        return redirect(url_for("main.index"))
    # get workbook names from list
    workbooks_names = [wb.name for wb in verification.data_export_request.workbooks]
    workbooks = ", ".join(workbooks_names)
    sas_url = services.data_export.export.make_sas_link(
        verification.data_export_request
    )

    # record export
    message = services.data_export_history.DataExportMessage(
        verification.data_export_request.requestor_person.user.id,
        verification.data_export_request.WorkGroup.id,
        # take the first workbook as exports from multiple is not supported
        [wb.id for wb in verification.data_export_request.workbooks][0],
        [reaction.id for reaction in verification.data_export_request.reactions],
        datetime.now().strftime("%Y-%m-%d"),
    )
    services.data_export_history.send_message(message)

    return render_template(
        "data_export/data_export_download.html",
        sas_url=sas_url,
        data_export_request=verification.data_export_request,
        workbooks=workbooks,
        workgroups=get_workgroups(),
        notification_number=get_notification_number(),
    )


@export_data_bp.route("/export_permission", methods=["GET", "POST"])
@login_required
@export_data_bp.doc(security="sessionAuth")
def export_permission():
    """
    Checks if the user has permission to export from the selected workbook.
    The requestor must be either a workbook member or a workgroup PI

    Returns:
        flask.Response: A JSON response with the status of the request
    """
    if security_pi_workgroup(
        request.json.get("workgroup") and request.json.get("workbook")
    ) or security_member_workgroup_workbook(
        request.json.get("workgroup"), request.json.get("workbook")
    ):
        return jsonify("permission accepted")
    else:
        return jsonify("permission denied")


@export_data_bp.route("/get_reaction_id_list", methods=["GET", "POST"])
@login_required
@export_data_bp.doc(security="sessionAuth")
def get_reaction_id_list():
    """
    Get a list of reaction IDs from a workbook

    Returns:
        flask.Response: A JSON response with the list of reaction
    """
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        request.json["workgroup"], request.json["workbook"]
    )
    if workbook:
        validated_reactions = [
            rxn
            for rxn in workbook.reactions
            if services.data_export.utils.validate_reaction(rxn)
        ]
        reaction_ids = [rxn.reaction_id for rxn in validated_reactions]
    else:
        reaction_ids = []
    return jsonify(reaction_ids)
