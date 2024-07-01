from datetime import datetime

from flask import redirect, request # renders html templates
from flask import Response, flash, render_template, url_for
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from sources.decorators import principal_investigator_required
from sources import models, services
from sources.auxiliary import get_notification_number, get_workgroups
from sources.extensions import db

from . import join_workgroup_bp  # imports the blueprint of the dummy route


@join_workgroup_bp.route("/join_workgroup/<workgroup>", methods=["GET", "POST"])
@login_required
def join_workgroup(workgroup=None) -> Response:
    # must be logged in
    workgroups_dropdown = get_workgroups()
    notification_number = get_notification_number()

    if workgroup in workgroups_dropdown:
        flash("You are already a member of this workgroup!")
        return redirect(url_for("main.index"))
    # find all the PIs for the workgroup
    pis = (
        db.session.query(models.Person)
        .join(models.Person.workgroup_principal_investigator)
        .filter(models.WorkGroup.name == workgroup)
        .all()
    )

    wg = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )

    # find person
    person = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .first()
    )

    # check to see if there is already an active request to join this workgroup from this user
    duplicates = (
        db.session.query(models.WGStatusRequest)
        .join(models.Person, models.WGStatusRequest.person == models.Person.id)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.id == wg.id)
        .filter(models.WGStatusRequest.status == 'active')
        .all()
    )
    if duplicates:
        flash(
            "You have already submitted a membership request for this workgroup. You will receive a notification "
            "when your request has been considered."
        )
        return redirect(url_for("main.index"))
    for pi in pis:
        notification = models.Notification(
            person=pi.id,
            type="New Workgroup Membership Request",
            info="You have a new request for a member to join Workgroup, "
            + workgroup
            + ", of which you are Principal Investigator.",
            time=datetime.now(),
            status="active",
            wg=wg.name,
            wb="",
            wg_request="",
        )
        services.email.send_notification(pi)
        db.session.add(notification)
        db.session.commit()  # Have to commit to get the id for notification
        wg_request = models.WGStatusRequest(
            principal_investigator=pi.id,
            person=person.id,
            wg=wg.id,
            current_role="Non-Member",
            new_role="Standard Member",
            time=datetime.now(),
            status="active",
            notification=notification.id,
        )
        db.session.add(wg_request)
    db.session.commit()
    flash(
        "Your membership has been requested. You will receive a notification when your request has been considered."
    )
    return redirect(url_for("main.index"))
# return render_template(
    #     "join_workgroup.html",
    #     form=form,
    #     workgroups=workgroups_dropdown,
    #     notification_number=notification_number,
    # )
