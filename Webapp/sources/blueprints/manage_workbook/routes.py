from datetime import datetime
from typing import Optional as OptionalType

from flask import Response, flash, jsonify, redirect, render_template, url_for
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from flask_wtf import FlaskForm
from sources import models, services
from sources.auxiliary import (
    duplicate_notification_check,
    get_all_workgroup_members,
    get_notification_number,
    get_workgroups,
    security_member_workgroup,
    security_pi_sr_workgroup,
)
from sources.extensions import db
from wtforms import SelectField, SubmitField
from wtforms.validators import Optional

from . import manage_workbook_bp  # imports the blueprint of the dummy route


class SelectWorkbookForm(FlaskForm):
    workbooks = SelectField("Workbook", coerce=int, validators=[Optional()])
    workbook_users = SelectField("Workbook members", validators=[Optional()])
    submit_wb = SubmitField(
        "Remove user from workbook", render_kw={"class": "btn btn-primary"}
    )


class JoinWorkbookForm(FlaskForm):
    workbooks = SelectField("Workbook", coerce=str)
    submit = SubmitField(
        "Request to Join Workbook", render_kw={"class": "btn btn-primary"}
    )


# Access workbook management page
@manage_workbook_bp.route("/manage_workbook/<workgroup>", methods=["GET", "POST"])
@manage_workbook_bp.route(
    "/manage_workbook/<workgroup>/<has_request>", methods=["GET", "POST"]
)
@manage_workbook_bp.route(
    "/manage_workbook/<workgroup>/<has_request>/<workbook>", methods=["GET", "POST"]
)
@login_required
def manage_workbook(
    workgroup: str, has_request: str = "no", workbook: OptionalType[str] = None
) -> Response:
    # must be logged in and a PI or SR of the workgroup
    if not security_pi_sr_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    current_workgroup = workgroup
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    # check if user is PI or SR and able to access this page
    pi = (
        db.session.query(models.User)
        .join(models.Person)
        .join(models.t_Person_WorkGroup)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == current_workgroup)
        .all()
    )
    pi_check = [x.email for x in pi]
    sr = (
        db.session.query(models.User)
        .join(models.Person)
        .join(models.t_Person_WorkGroup_2)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == current_workgroup)
        .all()
    )
    sr_check = [x.email for x in sr]

    if current_user.email not in pi_check and current_user.email not in sr_check:
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    """This function provides the initial dropdown choices for a user and then handles submission of the form"""
    # Initiate form
    form = SelectWorkbookForm()

    # Find all workbooks a user is a PI in
    wb1 = (
        db.session.query(models.WorkBook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == current_workgroup)
        .join(models.t_Person_WorkGroup)
        .join(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .all()
    )
    workbooks1 = [(wb.id, wb.name) for wb in wb1]

    # Find all workbooks a user is a senior researcher in
    wb2 = (
        db.session.query(models.WorkBook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == current_workgroup)
        .join(models.t_Person_WorkGroup_2)
        .join(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .all()
    )
    workbooks2 = [(wb.id, wb.name) for wb in wb2]
    form.workbooks.choices = workbooks1 + workbooks2
    if not workbook:
        try:
            workbook = form.workbooks.choices[0][0]
        except IndexError:
            return render_template(
                "manage_workbook.html",
                form=form,
                workgroups=workgroups,
                workbook=workbook,
                current_workgroup=current_workgroup,
                has_request="no",
                notification_number=notification_number,
            )
    form.workbooks.default = workbook
    form.process()

    # get all users and not users of workbook
    members = (
        db.session.query(models.Person)
        .join(models.t_Person_WorkBook)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook)
        .all()
    )
    wg = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == current_workgroup)
        .first()
    )
    pi, sr, sm = get_all_workgroup_members(wg)
    all_members = pi + sr + sm
    other_members = [x for x in all_members if x not in members]

    # get requests
    current_person = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .first()
    )
    wb = db.session.query(models.WorkBook).get(workbook)

    requests = (
        db.session.query(models.WBStatusRequest)
        .filter(models.WBStatusRequest.wb == wb.id)
        .filter(models.WBStatusRequest.status == "active")
        .filter(models.WBStatusRequest.pi_sr == current_person.id)
        .all()
    )
    return render_template(
        "manage_workbook.html",
        form=form,
        workgroups=workgroups,
        workbook=workbook,
        current_workgroup=current_workgroup,
        has_request=has_request,
        members=members,
        other_members=other_members,
        notification_number=notification_number,
        requests=requests,
    )


@manage_workbook_bp.route(
    "/manage_workbook/<workgroup>/<workbook>/<email>/<mode>", methods=["GET", "POST"]
)
@login_required
def add_remove_user_from_workbook(
    workgroup: str, workbook: str, email: str, mode: str
) -> Response:
    # must be logged in and a PI or SR of the workgroup
    if not security_pi_sr_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    user = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == email)
        .first()
    )
    wb = db.session.query(models.WorkBook).get(workbook)
    if mode == "remove":
        # remove user
        user.workbook_user.remove(wb)
        db.session.commit()
        return jsonify({"feedback": "This user has been removed from the workbook!"})
    else:
        # add user
        user.workbook_user.append(wb)
        db.session.commit()
        return jsonify({"feedback": "This user has been added to this workbook!"})


@manage_workbook_bp.route(
    "/manage_workbook_request/<workgroup>/<workbook>/<email>/<mode>",
    methods=["GET", "POST"],
)
@login_required
def manage_workbook_request(
    workgroup: str, workbook: str, email: str, mode: str
) -> Response:
    # must be logged in and a PI or SR of the workgroup
    if not security_pi_sr_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    wb = db.session.query(models.WorkBook).get(workbook)

    # change request to inactive
    request_objs = (
        db.session.query(models.WBStatusRequest)
        .filter(models.WBStatusRequest.wb == wb.id)
        .join(models.Person, models.WBStatusRequest.person == models.Person.id)
        .join(models.User)
        .filter(models.User.email == email)
        .all()
    )
    for obj in request_objs:
        obj.status = "inactive"

    person = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == email)
        .first()
    )

    # change initial notification to inactive
    for obj in [x.Notification for x in request_objs]:
        obj.status = "inactive"

    if mode == "deny":
        # notification to requester
        notification = models.Notification(
            person=person.id,
            type=f"Your Request to join {wb.name}",
            info=f"Your request to join Workbook, {wb.name}, has been denied.",
            time=datetime.now(),
            status="active",
        )
        db.session.add(notification)
        db.session.commit()
        services.email.send_notification(person)
        # feedback return
        return jsonify({"feedback": "This user has not been added to this workbook!"})
    else:
        # add user
        wb = db.session.query(models.WorkBook).get(workbook)
        person.workbook_user.append(wb)

        # notification to requester
        notification = models.Notification(
            person=person.id,
            type=f"Your Request to join {wb.name}",
            info=f"Your request to join Workbook, {wb.name}, has been approved.",
            time=datetime.now(),
            status="active",
        )
        db.session.add(notification)
        db.session.commit()
        services.email.send_notification(person)
        # feedback return
        return jsonify({"feedback": "This user has been added to this workbook!"})


@manage_workbook_bp.route("/join_workbook/<workgroup>", methods=["GET", "POST"])
@login_required
def join_workbook(workgroup: str) -> Response:
    # must be logged in and a member of the workgroup
    if not security_member_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    workgroups = get_workgroups()
    notification_number = get_notification_number()

    # create the form
    form = JoinWorkbookForm()

    # find all the workbooks in the WorkGroup database
    workbooks = (
        db.session.query(models.WorkBook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .all()
    )
    form.workbooks.choices = [x.name for x in workbooks]

    # when it is validly submitted
    if form.validate_on_submit():
        # get workgroup selected
        workbook = form.workbooks.data
        # check if user is already part of workbook, if so flash message rerender template
        wb = (
            db.session.query(models.WorkBook)
            .filter(models.WorkBook.name == workbook)
            .all()
        )
        email_list = [user.user.email for workbook in wb for user in workbook.users]

        if current_user.email in email_list:
            flash("You are already a member of this workbook!")
            return render_template(
                "join_workbook.html",
                form=form,
                workgroups=workgroups,
                notification_number=notification_number,
            )

        # find all the PI/SR for the workgroup
        pis = (
            db.session.query(models.Person)
            .join(models.Person.workgroup_principal_investigator)
            .filter(models.WorkGroup.name == workgroup)
            .all()
        )
        srs = (
            db.session.query(models.Person)
            .join(models.Person.workgroup_senior_researcher)
            .filter(models.WorkGroup.name == workgroup)
            .all()
        )
        pi_sr = pis + srs
        wb = (
            db.session.query(models.WorkBook)
            .filter(models.WorkBook.name == workbook)
            .join(models.WorkGroup)
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
        # set up notification and request for each one
        if duplicate_notification_check(
            [person],
            "New Workbook Membership Request",
            "active",
            workgroup,
            WB=workbook,
        ):
            flash(
                "You have already submitted a membership request for this workbook. You will receive a notification "
                "when your request has been considered."
            )
        for p in pi_sr:
            notification = models.Notification(
                person=p.id,
                type="New Workbook Membership Request",
                info=f"You have a new request for a member to join Workbook, {workbook}, of which you are Senior Researcher or Principal Investigator.",
                time=datetime.now(),
                status="active",
                wb=workbook,
                wg=workgroup,
            )
            services.email.send_notification(p)
            db.session.add(notification)
            db.session.commit()
            wb_request = models.WBStatusRequest(
                pi_sr=p.id,
                person=person.id,
                wb=wb.id,
                current_role="Non-Member",
                new_role="Member",
                time=datetime.now(),
                status="active",
                notification=notification.id,
            )
            db.session.add(wb_request)
        db.session.commit()
        flash(
            "Your membership has been requested. You will receive a notification when your request has been considered."
        )
        return redirect(url_for("main.index"))
    return render_template(
        "join_workbook.html",
        form=form,
        workgroups=workgroups,
        notification_number=notification_number,
    )


# from notification go to requests
@manage_workbook_bp.route(
    "/manage_workbook/go_to_workbook/<workbook>/<workgroup>", methods=["GET", "POST"]
)
@login_required
def go_to_workgroup(workbook: str, workgroup: str) -> Response:
    wb = (
        db.session.query(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )
    return redirect(
        url_for(
            "manage_workbook.manage_workbook",
            has_request="yes",
            workbook=wb.id,
            workgroup=workgroup,
        )
    )
