import threading

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
from flask_login import current_user, login_required
from sources import models, services
from sources.auxiliary import (
    get_notification_number,
    get_workgroups,
    security_member_workgroup_workbook,
    security_pi_workgroup,
)

from . import export_data_bp


# ### structure - 1 route per format
@export_data_bp.route("/export_data/home", methods=["GET", "POST"])
@login_required
def export_data_home():
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "data_export/data_export_home.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@export_data_bp.route("/export_data/new_request", methods=["POST", "GET"])
@login_required
def new_data_export_request():
    """Initiates a data export request if user has permission"""
    export_request = services.data_export.requests.NewRequest()
    permission_status = export_request.check_permissions()
    if permission_status == "permission accepted":
        export_request.submit_request()
    return jsonify(permission_status)


# route to reset password when link in email is used
@export_data_bp.route("/export_data/request_response/<token>", methods=["GET", "POST"])
@login_required
def export_data_request_response(token: str) -> Response:
    """
    Verify the token and the approver and either redirect or render the page with the request information
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
def export_denied():
    """Update the request in the database"""
    data_export_request = models.DataExportRequest.query.get(request.json["exportID"])
    request_status = services.data_export.requests.RequestStatus(data_export_request)
    request_status.deny()
    flash("Data export denied!")
    return redirect(url_for("main.index"))


@export_data_bp.route("/export_data/export_approved", methods=["POST", "GET"])
@login_required
def export_approved():
    """Update the database with the approval. Then check if all approvers have done so and start export."""
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
def request_download(token: str) -> Response:
    """
    Verify the token and the requestor and either redirect or render the page with the link to download the export file.
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
def export_permission():
    """
    Checks if the user has permission to export from the selected workbook.
    The requestor must be either a workbook member or a workgroup PI
    """
    if security_pi_workgroup(
        request.json["workgroup"]
    ) or security_member_workgroup_workbook(
        request.form["workgroup"], request.form["workbook"]
    ):
        return jsonify("permission accepted")
    else:
        return jsonify("permission denied")


@export_data_bp.route("/get_reaction_id_list", methods=["GET", "POST"])
@login_required
def get_reaction_id_list():
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        request.json["workgroup"], request.json["workbook"]
    )
    validated_reactions = [
        rxn
        for rxn in workbook.reactions
        if services.data_export.utils.validate_reaction(rxn)
    ]
    reaction_ids = [rxn.reaction_id for rxn in validated_reactions]
    return jsonify(reaction_ids)


#
# @export_data_bp.route("/export_data_eln_file", methods=["GET", "POST"])
# @login_required
# def export_data_eln_file() -> Response:
#     """
#     We are making a .eln file for a specific workbook. This is a zipped directory containing a ro-crate-metadata.json
#     to describe the contents.
#     """
#
#     # reaction_list = services.reaction.list_active_in_workbook(
#     #     workbook, workgroup, sort_crit="time"
#     # )
#     workgroup_name = request.form["workgroup"]
#     workbook_name = request.form["workbook"]
#     eln_file = ELNFile(workgroup_name, workbook_name)
#     # eln_file = DataExport(workgroup_name, workbook_name).to_eln_file()
#     return eln_file

# made according to this specification https://github.com/TheELNConsortium/TheELNFileFormat
# first we describe the ro-crate metadata_json
# ro_crate_metadata_json_contents = describe_ro_crate_metadata_json()
# now for each reaction we want to make a research object crate
# for idx, reaction in enumerate(reaction_list):
# make a folder per experiment.
