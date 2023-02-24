from flask import render_template, flash, redirect, url_for  # renders html templates
from flask_login import login_required  # protects a view function against anonymous users
from pony.orm import select
from flask_login import current_user
from sources.join_workgroup import join_workgroup_bp  # imports the blueprint of the dummy route
from sources import db
from sources.join_workgroup.forms import JoinWorkgroupForm
from sources.auxiliary import get_workgroups, get_notification_number, duplicate_notification_check
from datetime import datetime
from sources.email_methods import send_notification_email


# put routes here
@join_workgroup_bp.route('/join_workgroup', methods=['GET', 'POST'])
@login_required
def join_workgroup():
    # must be logged in
    workgroups_dropdown = get_workgroups()
    notification_number = get_notification_number()
    # create the form
    form = JoinWorkgroupForm()
    # find all the workgroups in the WorkGroup database
    workgroups = select(u for u in db.WorkGroup if u.principal_investigator)
    # add to list and add choices to form
    workgroup_list = []
    for name in workgroups:
        workgroup_list.append(name.name)
    form.workgroups.choices = workgroup_list
    # when it is validly submitted
    if form.validate_on_submit():
        # get workgroup selected
        workgroup = form.workgroups.data
        # check if user is already part of workgroup, if so flash message rerender template
        if workgroup in workgroups_dropdown:
            flash("You are already a member of this workgroup!")
            return render_template('join_workgroup.html', form=form, workgroups=workgroups,
                                   notification_number=notification_number)
        # find all the PIs for the workgroup
        pis = select(x.principal_investigator for x in db.WorkGroup if x.name == workgroup)[:]
        wg = select(x for x in db.WorkGroup if x.name == workgroup).first()
        # find person
        person = select(p for p in db.Person if p.user.email == current_user.email).first()
        # check to see if there are any new requests
        if duplicate_notification_check([person], "New Workgroup Membership Request", "active", wg.name):
            flash('You have already submitted a membership request for this workgroup. You will receive a notification '
                  'when your request has been considered.')
            return redirect(url_for('main.index'))
        for pi in pis:
            notification = db.Notification(person=pi, type="New Workgroup Membership Request",
                                           info="You have a new request for a member to join Workgroup, "
                                                + workgroup + ", of which you are Principal Investigator.",
                                           time=datetime.now(), status="active", WG=wg.name)
            send_notification_email(pi)
            db.WGStatusRequest(principal_investigator=pi, person=person, WG=wg, current_role="Non-Member",
                               new_role="Standard Member", time=datetime.now(), status="active",
                               notification=notification)
        flash('Your membership has been requested. You will receive a notification when your request has been considered.')
        return redirect(url_for('main.index'))
    return render_template('join_workgroup.html', form=form, workgroups=workgroups_dropdown,
                           notification_number=notification_number)



