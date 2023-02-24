from flask import render_template, jsonify, flash, redirect, url_for
from flask_login import login_required, current_user  # protects a view function against anonymous users
from sources.manage_workgroup import manage_workgroup_bp
from pony.orm import select
from sources import db
from sources.auxiliary import get_workgroups, get_notification_number, security_pi_workgroup, \
    security_member_workgroup, duplicate_notification_check,make_objects_inactive, workgroup_dict, get_all_workgroup_member_types
from datetime import datetime
from sources.email_methods import send_notification_email


@manage_workgroup_bp.route('/manage_workgroup/<workgroup>', methods=['GET', 'POST'])
@manage_workgroup_bp.route('/manage_workgroup/<workgroup>/<has_request>', methods=['GET', 'POST'])
@login_required
def manage_workgroup(workgroup, has_request="no"):
    # must be logged in and a PI of the workgroup
    if not security_pi_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    workgroups = get_workgroups()
    current_workgroup = workgroup
    # current_workgroup = str(request.form['workgroup'])
    notification_number = get_notification_number()
    # get all the pis, srs, and sms in the group
    wg = select(x for x in db.WorkGroup if x.name == current_workgroup).first()
    pi, sr, sm = get_all_workgroup_member_types(wg)
    # get requests for this workgroup
    current_person = select(x for x in wg.principal_investigator if x.user.email == current_user.email).first()
    requests = select(x for x in db.WGStatusRequest if x.WG == wg and x.principal_investigator == current_person and
                      x.status == "active")[:]
    return render_template('manage_workgroup.html', workgroups=workgroups, current_workgroup=current_workgroup, pi=pi,
                           sr=sr, sm=sm, requests=requests, has_request=has_request,
                           notification_number=notification_number)


@manage_workgroup_bp.route('/manage_workgroup/make_change/<workgroup>/<email>/<mode>/<current_status>',  methods=['GET', 'POST'])
@login_required
def make_change_to_workgroup(workgroup, email, mode, current_status):
    # must be logged in and a PI of the workgroup
    if not security_pi_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    # find the person and workgroup objects that the change relates to
    wg = select(x for x in db.WorkGroup if x.name == workgroup).first()
    person = select(x for x in db.Person if x.user.email == email).first()
    # mode is either remove user from workgroup or change role in workgroup
    if mode == "remove":
        if current_status == 'pi' and len(db.WorkGroup[wg.id].principal_investigator) < 2:
            return jsonify({"feedback": "Cannot remove PI. Workgroups must have at least one PI."})
        else:
            # uses the current_status variable as a dict key to determine the attribute being removed.
            getattr(db.Person[person.id], workgroup_dict[current_status]["person_to_wg_attr"]).remove(db.WorkGroup[wg.id])
            return jsonify({"feedback": f"{person.user.fullname} "
                                        f"has been removed from {wg.name} and associated Workbooks"})
    # else mode is changing role
    else:
        old_role, new_role = mode.split("-to-")
        if old_role == "pi" and len(db.WorkGroup[wg.id].principal_investigator) < 2:
            return jsonify({"feedback": "Cannot remove PI. Workgroups must have at least one PI."})
        # remove as old role then add as new one
        getattr(db.Person[person.id], workgroup_dict[old_role]['person_to_wg_attr']).remove(db.WorkGroup[wg.id])
        getattr(db.Person[person.id], workgroup_dict[new_role]['person_to_wg_attr']).add(db.WorkGroup[wg.id])
        return jsonify({"feedback": f"{person.user.fullname} has changed role from "
                                    f"{workgroup_dict[old_role]['display_string']} to "
                                    f"{workgroup_dict[new_role]['display_string']}"})


# PI/SR status request page
@manage_workgroup_bp.route('/PI_status_request/<user_type>/<new_role>/<workgroup>', methods=['GET', 'POST'])
@manage_workgroup_bp.route('/SR_status_request/<user_type>/<new_role>/<workgroup>', methods=['GET', 'POST'])
@login_required
def status_request(user_type, new_role, workgroup):
    # must be logged in and a member of the workgroup
    if not security_member_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    # find all the PIs for the workgroup
    pis = select(x.principal_investigator for x in db.WorkGroup if x.name == workgroup)[:]
    wg = select(x for x in db.WorkGroup if x.name == workgroup).first()
    # find the person
    person = select(x for x in db.Person if x.user.email == current_user.email).first()
    if duplicate_notification_check([person], "New Workgroup Role Reassignment Request", "active", wg.name):
        flash('You have already submitted a request for this workgroup. You will receive a notification '
              'when your request has been considered.')
        return redirect(url_for('main.index'))
    # set up notification and request for each one
    for pi in pis:
        notification = db.Notification(person=pi, type="New Workgroup Role Reassignment Request",
                                       info="You have a new request from a member to change their role in Workgroup, "
                                            + workgroup + ", of which you are Principal Investigator.",
                                       time=datetime.now(), status="active", WG=wg.name)
        send_notification_email(pi)
        db.WGStatusRequest(principal_investigator=pi, person=person, WG=wg, current_role=user_type,
                           new_role=new_role, time=datetime.now(), status="active", notification=notification)
    flash("Your request has been made. You will receive a notification when your request has been considered.")
    return redirect(url_for("main.index"))


# PI/SR status request decision
@manage_workgroup_bp.route('/manage_workgroup/change_status_request/<current_workgroup>/<email>/<mode>/<decision>', methods=['GET', 'POST'])
@login_required
def change_status_from_request(current_workgroup, email, mode, decision):
    # must be logged in and a PI of the workgroup
    if not security_pi_workgroup(current_workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    wg = select(x for x in db.WorkGroup if x.name == current_workgroup).first()
    # change request to inactive
    request_objs = select(x for x in db.WGStatusRequest if x.WG == wg and x.person.user.email == email)[:]
    make_objects_inactive(request_objs)
    # change initial notification to inactive
    notification_objs = select(x.notification for x in db.WGStatusRequest if x.WG == wg and x.person.user.email == email)[:]
    make_objects_inactive(notification_objs)
    person = select(x for x in db.Person if x.user.email == email).first()

    if decision == "deny":
        # send new notification
        if mode == "to-sm":
            db.Notification(person=person, type="Your Request to join " + wg.name,
                            info="Your request to join Workgroup, " + wg.name +
                                 ", has been denied.", time=datetime.now(), status="active")
            send_notification_email(person)
        elif mode == "sr-to-pi" or "sm-to-pi":
            db.Notification(person=person, type="Your Request to become a Principal Investigator",
                            info="Your request to become a Principal Investigator in Workgroup, " + wg.name +
                                 ", has been denied.", time=datetime.now(), status="active")
            send_notification_email(person)
        else:
            db.Notification(person=person, type="Your Request to become a Senior Researcher",
                            info="You request to become a Senior Researcher in Workgroup, " + wg.name +
                                 ", has been denied.", time=datetime.now(), status="active")
            send_notification_email(person)
        return jsonify({"feedback": "This request has been denied!"})
    # decision is to approve user
    else:
        if len(mode.split("-to-")) > 1:
            old_role, new_role = mode.split("-to-")
            # remove old role
            getattr(db.Person[person.id], workgroup_dict[old_role]['person_to_wg_attr']).remove(db.WorkGroup[wg.id])
            # set variables according to new role and workgroup
            new_role_display_string = workgroup_dict[new_role]["display_string"]
            request_type = f"Request to become a {new_role_display_string}"
            info = f"Your request to become a {new_role_display_string} in Workgroup, {wg.name}, " \
                   f"has been approved!"
            feedback = f"{person.user.fullname} has changed role from {workgroup_dict[old_role]['display_string']} " \
                       f"to {workgroup_dict[new_role]['display_string']}"
        else:
            # user is joining workgroup and has no old role
            new_role = "sm"
            request_type = "Request to join a Workgroup"
            info = f"Your request to join Workgroup {wg.name}, has been approved!"
            feedback = f"{person.user.fullname} has been added to {wg.name} as a Standard Member!"
        # add person to new role and send notification
        getattr(db.Person[person.id], workgroup_dict[new_role]['person_to_wg_attr']).add(db.WorkGroup[wg.id])
        db.Notification(person=person, type=request_type, info=info, time=datetime.now(), status="active")
        send_notification_email(person)
        return jsonify({"feedback": feedback})


# from notification go to requests
@manage_workgroup_bp.route('/manage_workgroup/go_to_workgroup/<workgroup>', methods=['GET', 'POST'])
@login_required
def go_to_workgroup(workgroup):
    return redirect(url_for("manage_workgroup.manage_workgroup", has_request="yes", workgroup=workgroup))
