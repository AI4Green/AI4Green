from flask import render_template, jsonify, redirect, url_for, flash
from flask_wtf import FlaskForm
from wtforms import SubmitField, SelectField
from wtforms.validators import Optional
from flask_login import login_required, current_user  # protects a view function against anonymous users
from sources.manage_workbook import manage_workbook_bp  # imports the blueprint of the dummy route
from pony.orm import select
from datetime import datetime
from sources import db
from sources.auxiliary import get_workgroups, get_notification_number, security_pi_sr_workgroup, \
    security_member_workgroup, duplicate_notification_check, get_all_workgroup_member_types
from sources.email_methods import send_notification_email


class SelectWorkbookForm(FlaskForm):
    workbooks = SelectField('Workbook', coerce=int, validators=[Optional()])
    workbook_users = SelectField('Workbook members', validators=[Optional()])
    submit_wb = SubmitField('Remove user from workbook', render_kw={"class": "btn btn-primary"})


class JoinWorkbookForm(FlaskForm):
    workbooks = SelectField('Workbook', coerce=str)
    submit = SubmitField('Request to Join Workbook', render_kw={"class": "btn btn-primary"})


# Access workbook management page
@manage_workbook_bp.route('/manage_workbook/<workgroup>', methods=['GET', 'POST'])
@manage_workbook_bp.route('/manage_workbook/<workgroup>/<has_request>', methods=['GET', 'POST'])
@manage_workbook_bp.route('/manage_workbook/<workgroup>/<has_request>/<workbook>', methods=['GET', 'POST'])
@login_required
def manage_workbook(workgroup, has_request="no", workbook=None):
    # must be logged in and a PI or SR of the workgroup
    if not security_pi_sr_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    current_workgroup = workgroup
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    # check if user is PI or SR and able to access this page
    pi_check = select(x.principal_investigator.user.email for x in db.WorkGroup if x.name == current_workgroup)[:]
    sr_check = select(x.senior_researcher.user.email for x in db.WorkGroup if x.name == current_workgroup)[:]
    if current_user.email not in pi_check and current_user.email not in sr_check:
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    """This function provides the initial dropdown choices for a user and then handles submission of the form"""
    # Initiate form
    form = SelectWorkbookForm()
    # Find all workbooks a user is a PI in
    workbooks1 = list(select((wb.id, wb.name) for wb in db.WorkBook if current_user.email in
                             wb.group.principal_investigator.user.email and wb.group.name == current_workgroup)[:])
    # Find all workbooks a user is a senior researcher in
    workbooks2 = list(select((wb.id, wb.name) for wb in db.WorkBook if current_user.email in
                             wb.group.senior_researcher.user.email and wb.group.name == current_workgroup)[:])
    form.workbooks.choices = workbooks1 + workbooks2
    if not workbook:
        try:
            workbook = form.workbooks.choices[0][0]
        except IndexError:
            return render_template('manage_workbook.html', form=form, workgroups=workgroups, workbook=workbook,
                                   current_workgroup=current_workgroup, has_request="no",
                                   notification_number=notification_number)
    form.workbooks.default = workbook
    form.process()
    # get all users and not users of workbook
    members = select(p for p in db.Person if workbook in p.workbook_user.id)[:]
    wg = select(x for x in db.WorkGroup if x.name == current_workgroup).first()
    pi, sr, sm = get_all_workgroup_member_types(wg)
    all_members = pi + sr + sm
    other_members = [x for x in all_members if x not in members]
    # get requests
    current_person = select(x for x in db.Person if x.user.email == current_user.email).first()
    wb = db.WorkBook[workbook]
    requests = select(x for x in db.WBStatusRequest if x.WB == wb and x.pi_sr == current_person and
                      x.status == "active")[:]
    return render_template('manage_workbook.html', form=form, workgroups=workgroups, workbook=workbook,
                           current_workgroup=current_workgroup, has_request=has_request,
                           members=members, other_members=other_members, notification_number=notification_number,
                           requests=requests)


@manage_workbook_bp.route('/manage_workbook/<workgroup>/<workbook>/<email>/<mode>', methods=['GET', 'POST'])
@login_required
def add_remove_user_from_workbook(workgroup, workbook, email, mode):
    # must be logged in and a PI or SR of the workgroup
    if not security_pi_sr_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    if mode == "remove":
        user = select(p for p in db.Person if workbook in p.workbook_user.id and p.user.email == email).first()
        # remove user
        db.Person[user.id].workbook_user.remove(db.WorkBook[workbook])
        # feedback
        return jsonify({"feedback": "This user has been removed from the workbook!"})
        # return
    else:
        current_workgroup = workgroup
        wg = select(x for x in db.WorkGroup if x.name == current_workgroup).first()
        user = select(x for x in wg.principal_investigator if x.user.email == email).first()
        if not user:
            user = select(x for x in wg.senior_researcher if x.user.email == email).first()
        if not user:
            user = select(x for x in wg.standard_member if x.user.email == email).first()
        # add user
        db.Person[user.id].workbook_user.add(db.WorkBook[workbook])
        # feedback return
        return jsonify({"feedback": "This user has been added to this workbook!"})


@manage_workbook_bp.route('/manage_workbook_request/<workgroup>/<workbook>/<email>/<mode>', methods=['GET', 'POST'])
@login_required
def manage_workbook_request(workgroup, workbook, email, mode):
    # must be logged in and a PI or SR of the workgroup
    if not security_pi_sr_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    wb = db.WorkBook[workbook]
    # change request to inactive
    request_objs = select(x for x in db.WBStatusRequest if x.WB == wb and x.person.user.email == email)[:]
    for obj in request_objs:
        obj.status = "inactive"
        person = obj.person
    # change initial notification to inactive
    notification_objs = select(x.notification for x in db.WBStatusRequest if x.WB == wb and x.person.user.email == email)[:]
    for obj in notification_objs:
        obj.status = "inactive"
    if mode == "deny":
        # notification to requester
        db.Notification(person=person, type="Your Request to join " + wb.name,
                        info="Your request to join Workbook, " + wb.name +
                             ", has been denied.", time=datetime.now(), status="active")
        send_notification_email(person)
        # feedback return
        return jsonify({"feedback": "This user has not been added to this workbook!"})
    else:
        # add user
        db.Person[person.id].workbook_user.add(db.WorkBook[workbook])
        # notification to requester
        db.Notification(person=person, type="Your Request to join " + wb.name,
                        info="Your request to join Workbook, " + wb.name +
                             ", has been approved.", time=datetime.now(), status="active")
        send_notification_email(person)
        # feedback return
        return jsonify({"feedback": "This user has been added to this workbook!"})


@manage_workbook_bp.route('/join_workbook/<workgroup>', methods=['GET', 'POST'])
@login_required
def join_workbook(workgroup):
    # must be logged in and a member of the workgroup
    if not security_member_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    # create the form
    form = JoinWorkbookForm()
    # find all the workbooks in the WorkGroup database
    workbooks = select(u.name for u in db.WorkBook if u.group.name == workgroup)[:]
    form.workbooks.choices = workbooks
    # when it is validly submitted
    if form.validate_on_submit():
        # get workgroup selected
        workbook = form.workbooks.data
        # check if user is already part of workbook, if so flash message rerender template
        email_list = select(u.users.user.email for u in db.WorkBook if u.name == workbook)[:]
        if current_user.email in email_list:
            flash("You are already a member of this workbook!")
            return render_template('join_workbook.html', form=form, workgroups=workgroups,
                                   notification_number=notification_number)
        # find all the PI/SR for the workgroup
        pis = select(x.principal_investigator for x in db.WorkGroup if x.name == workgroup)[:]
        srs = select(x.senior_researcher for x in db.WorkGroup if x.name == workgroup)[:]
        pi_sr = pis + srs
        wb = select(x for x in db.WorkBook if x.name == workbook and x.group.name == workgroup).first()
        # find person
        person = select(p for p in db.Person if p.user.email == current_user.email).first()
        # set up notification and request for each one
        if duplicate_notification_check([person], "New Workbook Membership Request", "active", workgroup, WB=workbook):
            flash('You have already submitted a membership request for this workbook. You will receive a notification '
                  'when your request has been considered.')
        for p in pi_sr:
            notification = db.Notification(person=p, type="New Workbook Membership Request",
                                           info="You have a new request for a member to join Workbook, "
                                                + workbook +
                                                ", of which you are Senior Researcher or Principal Investigator.",
                                           time=datetime.now(), status="active", WB=workbook, WG=workgroup)
            send_notification_email(p)
            db.WBStatusRequest(pi_sr=p, person=person, WB=wb, current_role="Non-Member", new_role="Member",
                               time=datetime.now(), status="active", notification=notification)
        flash('Your membership has been requested. You will receive a notification when your request has been considered.')
        return redirect(url_for('main.index'))
    return render_template('join_workbook.html', form=form, workgroups=workgroups,
                           notification_number=notification_number)


# from notification go to requests
@manage_workbook_bp.route('/manage_workbook/go_to_workbook/<workbook>/<workgroup>', methods=['GET', 'POST'])
@login_required
def go_to_workgroup(workbook, workgroup):
    wb = select(x for x in db.WorkBook if x.name == workbook and x.group.name == workgroup).first()
    return redirect(url_for("manage_workbook.manage_workbook", has_request="yes", workbook=wb.id, workgroup=workgroup))
