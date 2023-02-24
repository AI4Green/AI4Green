from flask import render_template, flash, redirect, url_for
from flask_login import login_required
from sources.admin_dashboard import admin_dashboard_bp
from sources import db
from datetime import datetime
from pony.orm import select
from sources.auxiliary import get_workgroups, get_notification_number, security_admin_only
from sources.email_methods import send_notification_email


@admin_dashboard_bp.route('/admin_dashboard', methods=['GET', 'POST'])
@admin_dashboard_bp.route('/admin_dashboard/<request_institution>/<request_name>/<decision>', methods=['GET', 'POST'])
@login_required
def admin_dashboard(request_institution=None, request_name=None, decision=None):
    """These routes display all active requests from database and perform the actions on approval/denial"""
    # must be logged in and admin
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    if not security_admin_only():
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    message = None
    # if block to handle request if present - otherwise just renders template with requests/reports
    if request_institution:
        pi_request = select(x for x in db.WorkGroup_request if x.institution.name == request_institution and x.name == request_name).first()
        approval_workgroup = pi_request.workgroup
        print(approval_workgroup.name)
        if pi_request.status == "inactive":
            message = "This request has already been considered"
        else:
            # set the request to inactive now it has been responded to
            pi_request.status = "inactive"
            if decision == "deny":
                db.Notification(person=pi_request.principal_investigator, type="Your Request to Create a Workgroup",
                                info="Your Workgroup request for " + pi_request.name + " at " +
                                     pi_request.institution.name + " has been denied.", time=datetime.now(),
                                status="active")
                send_notification_email(pi_request.principal_investigator)
                message = "This request has been denied and workgroup removed"
                approval_workgroup.delete()
            if decision == "approve":
                # check if the workgroup has already beeba approved and if so return msg advising admin to check the db
                workgroup_already_approved = select(u for u in db.WorkGroup if u.institution.name == pi_request.institution.name
                                                          and u.name == pi_request.name and u.approved is True).first()
                if workgroup_already_approved:
                    pi_request.status = "active"
                    message = "There was a problem with this request. Check the database for further details."
                # expected behaviour is to follow the else block as workgroup shouldn't be in requests if it is approved
                else:
                    # When request has been approved the PI receives a notification and the workgroup approval status is updated
                    db.Notification(person=pi_request.principal_investigator, type="Your Request to Create a Workgroup",
                                    info="Your Workgroup request for " + pi_request.name + " at " +
                                         pi_request.institution.name + " has been approved.", time=datetime.now(),
                                    status="active")
                    send_notification_email(pi_request.principal_investigator)
                    # update the workgroup approval status and flash message
                    approval_workgroup.approved = True
                    message = "This request has been approved"
    compound_data_error_reports = select(x for x in db.CompoundDataErrorReport)[:]  # exclude reports older than x months?
    number_compounds_in_db = len(select(x for x in db.Compound)[:])
    number_reactions_in_db = len(select(x for x in db.Reaction)[:])
    pi_requests = select(x for x in db.WorkGroup_request if x.status == "active")[:]
    users = select(x for x in db.User)[:]
    wgs = select(x for x in db.WorkGroup)[:]
    if message:
        flash(message)
    return render_template('admin_dashboard.html', workgroups=workgroups,
                           compound_data_error_reports=compound_data_error_reports,
                           number_compounds_in_db=number_compounds_in_db, number_reactions_in_db=number_reactions_in_db,
                           pi_requests=pi_requests, notification_number=notification_number, users=users, wgs=wgs)
