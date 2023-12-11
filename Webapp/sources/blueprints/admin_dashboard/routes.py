from datetime import datetime

from flask import Response, flash, redirect, render_template, url_for
from flask_login import login_required
from sources import models, services
from sources.auxiliary import (
    get_notification_number,
    get_workgroups,
    security_admin_only,
)
from sources.extensions import db

from . import admin_dashboard_bp


@admin_dashboard_bp.route("/admin_dashboard", methods=["GET", "POST"])
@admin_dashboard_bp.route(
    "/admin_dashboard/<request_institution>/<request_name>/<decision>",
    methods=["GET", "POST"],
)
@login_required
def admin_dashboard(
    request_institution: str = None, request_name: str = None, decision: str = None
) -> Response:
    """These routes display all active requests from database and perform the actions on approval/denial"""
    # must be logged in and admin
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    if not security_admin_only():
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))

    # if admin responds to a new workgroup request this block of code will process it and flash a message to the admin
    if request_institution:
        new_workgroup_request_message = handle_new_workgroup_request(
            request_institution, request_name, decision
        )
        flash(new_workgroup_request_message)

    # get data to populate the admin dashboard
    compound_data_error_reports = services.compound.get_compound_data_error_reports()
    reaction_count = services.reaction.count()
    compound_count = services.compound.count()
    new_workgroup_requests = services.workgroup.get_new_workgroup_requests()
    all_users = services.user.list_all()
    all_workgroups = services.workgroup.list_all()
    all_workbooks = services.workbook.list_all()
    recent_reactions = services.reaction.list_recent()

    return render_template(
        "admin_dashboard.html",
        workgroups=workgroups,
        compound_data_error_reports=compound_data_error_reports,
        number_compounds_in_db=compound_count,
        number_reactions_in_db=reaction_count,
        new_workgroup_requests=new_workgroup_requests,
        notification_number=notification_number,
        users=all_users,
        wgs=all_workgroups,
        workbooks=all_workbooks,
        recent_reactions=recent_reactions,
    )


def handle_new_workgroup_request(
    request_institution: str, request_name: str, decision: str
) -> str:
    """
    Processes the administrator's response to the new workgroup request. Either accepting or denying the request
    and updating the database accordingly.

    Args:
        request_institution - the name of the institution that the workgroup will be a part of
        request_name - the name of the new workgroup being requested
        decision - depends on button the admin clicks to either approve or deny the new workgroup

    Returns:
        message - to be flashed to admin to give feedback on the decision.

    """

    new_workgroup_request = (
        db.session.query(models.WorkGroupRequest)
        .join(models.WorkGroup)
        .filter(models.Institution.name == request_institution)
        .filter(models.WorkGroupRequest.name == request_name)
        .first()
    )
    if new_workgroup_request.status == "inactive":
        message = "This request has already been considered"
    else:
        message = "This request has not yet been considered"
        # set the request to inactive now it has been responded to
        new_workgroup_request.status = "inactive"
        approval_workgroup = new_workgroup_request.WorkGroup

        if decision == "deny":
            notification = models.Notification(
                person=new_workgroup_request.principal_investigator,
                type="Your Request to Create a Workgroup",
                info="Your Workgroup request for "
                + new_workgroup_request.name
                + " at "
                + new_workgroup_request.Institution.name
                + " has been denied.",
                time=datetime.now(),
                status="active",
                wg="",
                wb="",
                wg_request="",
            )
            db.session.add(notification)
            services.email.send_notification(new_workgroup_request.pi)
            message = "This request has been denied and workgroup removed"
            db.session.delete(approval_workgroup)
        if decision == "approve":
            # check if the workgroup has already been approved and if so return msg advising admin to check the db
            workgroup_already_approved = (
                db.session.query(models.WorkGroup)
                .filter(
                    models.Institution.name == new_workgroup_request.Institution.name
                )
                .filter(models.WorkGroup == new_workgroup_request.name)
                .filter(models.WorkGroup.approved is True)
                .first()
            )
            if workgroup_already_approved:
                new_workgroup_request.status = "active"
                message = "There was a problem with this request. Check the database for further details."
            # expected behaviour is to follow the else block as workgroup shouldn't be in requests if it is approved
            else:
                # When request has been approved the PI receives a notification and the workgroup approval status is updated
                notification = models.Notification(
                    person=new_workgroup_request.principal_investigator,
                    type="Your Request to Create a Workgroup",
                    info="Your Workgroup request for "
                    + new_workgroup_request.name
                    + " at "
                    + new_workgroup_request.Institution.name
                    + " has been approved.",
                    time=datetime.now(),
                    status="active",
                    wg="",
                    wb="",
                    wg_request="",
                )
                db.session.add(notification)
                services.email.send_notification(new_workgroup_request.pi)
                # update the workgroup approval status and flash message
                approval_workgroup.approved = True
                message = "This request has been approved"
        db.session.commit()
        return message
