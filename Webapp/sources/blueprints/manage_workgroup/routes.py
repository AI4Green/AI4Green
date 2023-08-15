from datetime import datetime

from flask import Response, flash, jsonify, redirect, render_template, url_for
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from sources import models, services
from sources.auxiliary import (
    duplicate_notification_check,
    get_all_workgroup_members,
    get_notification_number,
    get_workgroups,
    make_objects_inactive,
    security_member_workgroup,
    security_pi_workgroup,
    workgroup_dict,
)
from sources.extensions import db

from . import manage_workgroup_bp


@manage_workgroup_bp.route("/manage_workgroup/<workgroup>", methods=["GET", "POST"])
@manage_workgroup_bp.route(
    "/manage_workgroup/<workgroup>/<has_request>", methods=["GET", "POST"]
)
@login_required
def manage_workgroup(workgroup: str, has_request: str = "no") -> Response:
    # must be logged in and a PI of the workgroup
    if not security_pi_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    workgroups = get_workgroups()
    current_workgroup = workgroup
    notification_number = get_notification_number()
    # get all the pis, srs, and sms in the group
    wg = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == current_workgroup)
        .first()
    )
    pi, sr, sm = get_all_workgroup_members(wg)

    # get requests for this workgroup
    user = (
        db.session.query(models.Person)
        .join(models.WorkGroup.principal_investigator)
        .filter(
            models.WorkGroup.principal_investigator.any(
                models.Person.user.has(email=current_user.email)
            )
        )
        .first()
    )

    requests = (
        db.session.query(models.WGStatusRequest)
        .filter(models.WGStatusRequest.status == "active")
        .join(
            models.Person,
            models.WGStatusRequest.principal_investigator == models.Person.id,
        )
        .join(models.User)
        .filter(models.WGStatusRequest.principal_investigator == user.id)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.id == wg.id)
        .all()
    )
    return render_template(
        "manage_workgroup.html",
        workgroups=workgroups,
        current_workgroup=current_workgroup,
        pi=pi,
        sr=sr,
        sm=sm,
        requests=requests,
        has_request=has_request,
        notification_number=notification_number,
    )


@manage_workgroup_bp.route(
    "/manage_workgroup/make_change/<workgroup>/<email>/<mode>/<current_status>",
    methods=["GET", "POST"],
)
@login_required
def make_change_to_workgroup(
    workgroup: str, email: str, mode: str, current_status: str
) -> Response:
    # must be logged in and a PI of the workgroup
    if not security_pi_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    # find the person and workgroup objects that the change relates to
    wg = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )
    person = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == email)
        .first()
    )
    # mode is either remove user from workgroup or change role in workgroup
    if mode == "remove":
        if current_status == "pi" and len(wg.principal_investigator) < 2:
            return jsonify(
                {"feedback": "Cannot remove PI. Workgroups must have at least one PI."}
            )

        # uses the current_status variable as a dict key to determine the attribute being removed.
        getattr(person, workgroup_dict[current_status]["person_to_wg_attr"]).remove(wg)
        db.session.commit()
        return jsonify(
            {
                "feedback": f"{person.user.fullname} "
                f"has been removed from {wg.name} and associated Workbooks"
            }
        )
    else:
        old_role, new_role = mode.split("-to-")
        if old_role == "pi" and len(wg.principal_investigator) < 2:
            return jsonify(
                {"feedback": "Cannot remove PI. Workgroups must have at least one PI."}
            )
        # remove as old role then add as new one
        getattr(person, workgroup_dict[old_role]["person_to_wg_attr"]).remove(wg)
        getattr(person, workgroup_dict[new_role]["person_to_wg_attr"]).add(wg)
        db.session.commit()
        return jsonify(
            {
                "feedback": f"{person.user.fullname} has changed role from "
                f"{workgroup_dict[old_role]['display_string']} to "
                f"{workgroup_dict[new_role]['display_string']}"
            }
        )


# PI/SR status request page
@manage_workgroup_bp.route(
    "/PI_status_request/<user_type>/<new_role>/<workgroup>", methods=["GET", "POST"]
)
@manage_workgroup_bp.route(
    "/SR_status_request/<user_type>/<new_role>/<workgroup>", methods=["GET", "POST"]
)
@login_required
def status_request(user_type: str, new_role: str, workgroup: str) -> Response:
    # must be logged in and a member of the workgroup
    if not security_member_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    # find all the PIs for the workgroup
    pis = (
        db.session.query(models.Person)
        .join(models.t_Person_WorkGroup)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .all()
    )
    wg = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )

    # find the person
    person = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .first()
    )

    if duplicate_notification_check(
        [person], "New Workgroup Role Reassignment Request", "active", wg.name
    ):
        flash(
            "You have already submitted a request for this workgroup. You will receive a notification "
            "when your request has been considered."
        )
        return redirect(url_for("main.index"))

    # set up notification and request for each one
    for pi in pis:
        notification = models.Notification(
            person=pi.id,
            type="New Workgroup Role Reassignment Request",
            info=f"You have a new request from a member to change their role in Workgroup, {workgroup}, of which you are Principal Investigator.",
            time=datetime.now(),
            status="active",
            wg=wg.name,
        )
        db.session.add(notification)
        db.session.commit()  # needed to get the notification id.
        services.email.send_notification(pi)
        wg_request = models.WGStatusRequest(
            principal_investigator=pi.id,
            person=person.id,
            wg=wg.id,
            current_role=user_type,
            new_role=new_role,
            time=datetime.now(),
            status="active",
            notification=notification.id,
        )
        db.session.add(wg_request)
    db.session.commit()
    flash(
        "Your request has been made. You will receive a notification when your request has been considered."
    )
    return redirect(url_for("main.index"))


# PI/SR status request decision
@manage_workgroup_bp.route(
    "/manage_workgroup/change_status_request/<current_workgroup>/<email>/<mode>/<decision>",
    methods=["GET", "POST"],
)
@login_required
def change_status_from_request(
    current_workgroup: str, email: str, mode: str, decision: str
) -> Response:
    # must be logged in and a PI of the workgroup
    if not security_pi_workgroup(current_workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    wg = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == current_workgroup)
        .first()
    )

    # change request to inactive
    request_objs = (
        db.session.query(models.WGStatusRequest)
        .join(models.Person, models.WGStatusRequest.person == models.Person.id)
        .join(models.User)
        .filter(models.User.email == email)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.id == wg.id)
        .all()
    )
    make_objects_inactive(request_objs)

    # change initial notification to inactive
    notification_objs = (
        db.session.query(models.WGStatusRequest)
        .join(models.Person, models.WGStatusRequest.person == models.Person.id)
        .join(models.User)
        .filter(models.User.email == email)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.id == wg.id)
        .all()
    )
    notification_objs = [x.Notification for x in notification_objs]
    make_objects_inactive(notification_objs)

    person = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == email)
        .first()
    )

    if decision == "deny":
        # send new notification
        if mode == "to-sm":
            notification = models.Notification(
                person=person.id,
                type=f"Your Request to join {wg.name}",
                info=f"Your request to join Workgroup, {wg.name} has been denied.",
                time=datetime.now(),
                status="active",
            )
        elif mode == "sr-to-pi" or "sm-to-pi":
            notification = models.Notification(
                person=person.id,
                type="Your Request to become a Principal Investigator",
                info=f"Your request to become a Principal Investigator in Workgroup, {wg.name}, has been denied.",
                time=datetime.now(),
                status="active",
            )
        else:
            notification = models.Notification(
                person=person.id,
                type="Your Request to become a Senior Researcher",
                info=f"You request to become a Senior Researcher in Workgroup, {wg.name}, has been denied.",
                time=datetime.now(),
                status="active",
            )
        services.email.send_notification(person)
        db.session.add(notification)
        db.session.commit()
        return jsonify({"feedback": "This request has been denied!"})
    # decision is to approve user
    else:
        if len(mode.split("-to-")) > 1:
            old_role, new_role = mode.split("-to-")
            # remove old role
            getattr(person, workgroup_dict[old_role]["person_to_wg_attr"]).remove(wg)
            # set variables according to new role and workgroup
            new_role_display_string = workgroup_dict[new_role]["display_string"]
            request_type = f"Request to become a {new_role_display_string}"
            info = (
                f"Your request to become a {new_role_display_string} in Workgroup, {wg.name}, "
                f"has been approved!"
            )
            feedback = (
                f"{person.user.fullname} has changed role from {workgroup_dict[old_role]['display_string']} "
                f"to {workgroup_dict[new_role]['display_string']}"
            )
        else:
            # user is joining workgroup and has no old role
            new_role = "sm"
            request_type = "Request to join a Workgroup"
            info = f"Your request to join Workgroup {wg.name}, has been approved!"
            feedback = f"{person.user.fullname} has been added to {wg.name} as a Standard Member!"
        # add person to new role and send notification
        getattr(person, workgroup_dict[new_role]["person_to_wg_attr"]).add(wg)
        notification = models.Notification(
            person=person.id,
            type=request_type,
            info=info,
            time=datetime.now(),
            status="active",
        )
        db.session.add(notification)
        db.session.commit()
        services.email.send_notification(person)
        return jsonify({"feedback": feedback})


# from notification go to requests
@manage_workgroup_bp.route(
    "/manage_workgroup/go_to_workgroup/<workgroup>", methods=["GET", "POST"]
)
@login_required
def go_to_workgroup(workgroup: str) -> Response:
    return redirect(
        url_for(
            "manage_workgroup.manage_workgroup", has_request="yes", workgroup=workgroup
        )
    )
