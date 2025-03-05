from datetime import datetime

from flask import Response, flash, redirect, url_for
from flask_login import current_user, login_required
from sources import models, services
from sources.auxiliary import get_workgroups, make_objects_inactive, workgroup_dict
from sources.extensions import db

from . import join_workgroup_bp


@join_workgroup_bp.route("/join_workgroup/<workgroup>", methods=["GET", "POST"])
@login_required
@join_workgroup_bp.doc("sessionAuth")
def join_workgroup(workgroup=None) -> Response:
    """
    Sends a notification to the principal investigator of the workgroup to request membership.
    Args:
        workgroup: The workgroup to join

    Returns:
        flask.Response: A redirect to the home page with a message indicating whether the request was successful.

    """
    workgroups_dropdown = get_workgroups()

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
        .filter(models.WGStatusRequest.status == "active")
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


@join_workgroup_bp.route(
    "/join_workgroup_response/<notification>/<decision>", methods=["GET", "POST"]
)
def join_workgroup_response(notification=None, decision=None):
    """
    Accept or reject a request by a user to join the workgroup
    Args:
        notification: The notification being responded to
        decision: The decision to accept or reject the request

    Returns:
        flask.Response: A redirect to the home page with a message indicating whether the request was successful.
    """
    if notification:
        notification = (
            db.session.query(models.Notification)
            .filter(models.Notification.id == notification)
            .first()
        )

        make_objects_inactive([notification])

        workgroup = services.workgroup.from_name(notification.wg)

        wg_request = (
            db.session.query(models.WGStatusRequest)
            .join(models.Person, models.WGStatusRequest.person == models.Person.id)
            .filter(models.WGStatusRequest.notification == notification.id)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .join(models.WorkGroup)
            .filter(models.WorkGroup.id == workgroup.id)
            .first()
        )
        make_objects_inactive([wg_request])

        if decision == "accept":
            person = services.person.from_current_user_email()

            user_type = wg_request.new_role

            getattr(person, workgroup_dict[user_type]["person_to_wg_attr"]).add(
                workgroup
            )
            flash(
                f"You have been successfully added to the Workgroup: {workgroup.name}"
            )

        else:
            flash(f"You have not been added to the Workgroup: {workgroup.name}")

        db.session.commit()

        return redirect(url_for("main.index"))
