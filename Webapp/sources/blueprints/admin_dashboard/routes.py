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
    message = None
    # if block to handle request if present - otherwise just renders template with requests/reports
    if request_institution:
        pi_request = (
            db.session.query(models.WorkGroupRequest)
            .join(models.WorkGroup)
            .filter(models.Institution.name == request_institution)
            .filter(models.WorkGroupRequest.name == request_name)
            .first()
        )
        if pi_request.status == "inactive":
            message = "This request has already been considered"
        else:
            # set the request to inactive now it has been responded to
            pi_request.status = "inactive"
            approval_workgroup = pi_request.WorkGroup

            if decision == "deny":
                notification = models.Notification(
                    person=pi_request.principal_investigator,
                    type="Your Request to Create a Workgroup",
                    info="Your Workgroup request for "
                    + pi_request.name
                    + " at "
                    + pi_request.Institution.name
                    + " has been denied.",
                    time=datetime.now(),
                    status="active",
                    wg="",
                    wb="",
                    wg_request="",
                )
                db.session.add(notification)
                services.email.send_notification(pi_request.pi)
                message = "This request has been denied and workgroup removed"
                db.session.delete(approval_workgroup)
            if decision == "approve":
                # check if the workgroup has already been approved and if so return msg advising admin to check the db
                workgroup_already_approved = (
                    db.session.query(models.WorkGroup)
                    .filter(models.Institution.name == pi_request.Institution.name)
                    .filter(models.WorkGroup == pi_request.name)
                    .filter(models.WorkGroup.approved is True)
                    .first()
                )
                if workgroup_already_approved:
                    pi_request.status = "active"
                    message = "There was a problem with this request. Check the database for further details."
                # expected behaviour is to follow the else block as workgroup shouldn't be in requests if it is approved
                else:
                    # When request has been approved the PI receives a notification and the workgroup approval status is updated
                    notification = models.Notification(
                        person=pi_request.principal_investigator,
                        type="Your Request to Create a Workgroup",
                        info="Your Workgroup request for "
                        + pi_request.name
                        + " at "
                        + pi_request.Institution.name
                        + " has been approved.",
                        time=datetime.now(),
                        status="active",
                        wg="",
                        wb="",
                        wg_request="",
                    )
                    db.session.add(notification)
                    services.email.send_notification(pi_request.pi)
                    # update the workgroup approval status and flash message
                    approval_workgroup.approved = True
                    message = "This request has been approved"
            db.session.commit()

    compound_data_error_reports = db.session.query(
        models.CompoundDataErrorReport
    ).all()  # exclude reports older than x months?
    number_reactions_in_db = db.session.query(models.Reaction).count()
    number_compounds_in_db = db.session.query(models.Compound).count()
    pi_requests = (
        db.session.query(models.WorkGroupRequest)
        .filter(models.WorkGroupRequest.status == "active")
        .all()
    )
    users = db.session.query(models.User).all()
    wgs = db.session.query(models.WorkGroup).all()
    if message:
        flash(message)
    return render_template(
        "admin_dashboard.html",
        workgroups=workgroups,
        compound_data_error_reports=compound_data_error_reports,
        number_compounds_in_db=number_compounds_in_db,
        number_reactions_in_db=number_reactions_in_db,
        pi_requests=pi_requests,
        notification_number=notification_number,
        users=users,
        wgs=wgs,
    )
