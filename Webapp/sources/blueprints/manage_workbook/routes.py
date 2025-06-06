from datetime import datetime
from typing import Optional as OptionalType

from flask import Response, flash, jsonify, redirect, render_template, url_for
from flask_login import current_user, login_required
from flask_wtf import FlaskForm
from sources import models, services
from sources.auxiliary import (
    get_all_workgroup_members,
    get_notification_number,
    get_workgroups,
)
from sources.decorators import (
    principal_investigator_or_senior_researcher_required,
    workgroup_member_required,
)
from sources.extensions import db
from wtforms import SelectField, SubmitField
from wtforms.validators import Optional

from . import manage_workbook_bp


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
@principal_investigator_or_senior_researcher_required
@manage_workbook_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator or Senior Researcher role in the specified workgroup.",
)
def manage_workbook(
    workgroup: str, has_request: str = "no", workbook: OptionalType[str] = None
) -> Response:
    """
    Renders the manage workbook page
    Args:
        workgroup: the workgroup the workbook belongs to
        has_request: whether the workbook has a request to respond to
        workbook: the workbook being managed

    Returns:
        flask.Response: The rendered template for managing the selected workbook.
    """
    current_workgroup = workgroup
    workgroups = get_workgroups()
    notification_number = get_notification_number()

    # This function provides the initial dropdown choices for a user and then handles submission of the form
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
                "work_structures/manage_workbook.html",
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
        "work_structures/manage_workbook.html",
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
@principal_investigator_or_senior_researcher_required
@manage_workbook_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator or Senior Researcher role in the specified workgroup.",
)
def add_remove_user_from_workbook(
    workgroup: str, workbook: str, email: str, mode: str
) -> Response:
    """
    Add or remove a user from a workbook

    Args:
        workgroup: the workgroup the workbook belongs to
        workbook: the workbook to add or remove the user from
        email: the email of the user to add or remove
        mode: whether to add or remove the user

    Returns:
        flask.Response: A JSON response with feedback on the action taken
    """

    user = (
        db.session.query(models.Person)
        .join(models.User)
        .filter(models.User.email == email)
        .first()
    )
    wb = db.session.query(models.WorkBook).get(workbook)
    workgroup = db.session.query(models.WorkGroup).filter_by(name=workgroup).first()
    if mode == "remove":
        # remove user
        user.workbook_user.remove(wb)
        db.session.commit()
        # record access change
        message = services.data_access_history.DataAccessMessage(
            user.id,
            workgroup.id,
            "Access",
            "No Access",
            datetime.now().strftime("%Y-%m-%d"),
            wb.id,
        )
        services.data_access_history.send_message(message)
        return jsonify({"feedback": "This user has been removed from the workbook!"})
    else:
        # add user
        user.workbook_user.append(wb)
        db.session.commit()
        # record access change
        message = services.data_access_history.DataAccessMessage(
            user.id, workgroup.id, "No Access", "Access", wb.id
        )
        services.data_access_history.send_message(message)
        return jsonify({"feedback": "This user has been added to this workbook!"})


@manage_workbook_bp.route(
    "/manage_workbook_request/<workgroup>/<workbook>/<email>/<mode>",
    methods=["GET", "POST"],
)
@login_required
@principal_investigator_or_senior_researcher_required
@manage_workbook_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator or Senior Researcher role in the specified workgroup.",
)
def manage_workbook_request(
    workgroup: str, workbook: str, email: str, mode: str
) -> Response:
    """
    Manage a request to join a workbook

    Args:
        workgroup: the workgroup the workbook belongs to
        workbook: the workbook to accept or deny the request to join
        email: the email of the user requesting to join
        mode: whether to accept or deny the request

    Returns:
        flask.Response: A JSON response with feedback on the action taken
    """
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
@workgroup_member_required
@manage_workbook_bp.doc(
    security="sessionAuth",
    description="Requires login and membership in the specified workgroup.",
)
def join_workbook(workgroup: str) -> Response:
    # TODO - is this function deprecated?
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
                "work_structures/join_workbook.html",
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
        if services.notifications.duplicate_notification_check(
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
            return redirect(url_for("main.index"))
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
        "work_structures/join_workbook.html",
        form=form,
        workgroups=workgroups,
        notification_number=notification_number,
    )


@manage_workbook_bp.route(
    "/manage_workbook/go_to_workbook/<workbook>/<workgroup>", methods=["GET", "POST"]
)
@login_required
@manage_workbook_bp.doc(security="sessionAuth")
def go_to_workgroup(workbook: str, workgroup: str) -> Response:
    """
    Redirects to the workbook management page from a notification

    Args:
        workbook: the workbook to manage
        workgroup: the workgroup the workbook belongs to

    Returns:
        flask.Response: A redirect to the workbook management page
    """
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
