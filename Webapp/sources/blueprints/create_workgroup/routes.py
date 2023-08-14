from datetime import datetime

import pytz
from flask import redirect  # renders html templates
from flask import Response, flash, render_template, request, url_for
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from flask_wtf import FlaskForm
from sources import auxiliary, models, services
from sources.auxiliary import get_workgroups, remove_spaces_and_dashes
from sources.extensions import db
from sqlalchemy import func
from wtforms import SelectField, StringField, SubmitField, TextAreaField
from wtforms.validators import Length

from . import create_workgroup_bp  # imports the blueprint of the dummy route


class CreateWorkgroupForm(FlaskForm):
    # institution = SelectField('Institution')
    workgroup = StringField("Workgroup", validators=[Length(min=4, max=320)])
    info = TextAreaField(
        "Further Information",
        render_kw={
            "placeholder": "I am a Principal Investigator researching...",
            "rows": "4",
            "cols": "50",
        },
    )
    submit = SubmitField("Create Workgroup", render_kw={"class": "btn btn-primary"})


@create_workgroup_bp.route("/create_workgroup", methods=["GET", "POST"])
@login_required
def create_workgroup() -> Response:
    # must be logged in
    workgroups = get_workgroups()
    form = CreateWorkgroupForm()
    # institutions removed for now
    # get all institutions (should PIs be associated with one or more institutions?)
    # institutions = sq.session.query(models.Institution).all()
    # add to institution options
    # form.institution.choices = institutions
    # on valid submit
    if request.method == "POST" and form.validate_on_submit():
        # institutions removed for now
        # institution = form.institution.data
        workgroup_name = auxiliary.sanitise_user_input(form.workgroup.data)
        info = form.info.data

        # find PI object from current user email
        PI = (
            db.session.query(models.Person)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .first()
        )

        # validate no other workgroups awaiting approval
        other_workgroups_waiting_approval = (
            db.session.query(models.WorkGroup)
            .filter(models.WorkGroup.approved is False)
            .filter(
                models.WorkGroup.principal_investigator.any(models.Person.id == PI.id)
            )
            .all()
        )

        workgroup_name_delimiters_removed = remove_spaces_and_dashes(workgroup_name)
        if other_workgroups_waiting_approval:
            flash(
                "You cannot create another workgroup until your first workgroup has been approved"
            )
            return render_template(
                "create_workgroup.html", form=form, workgroups=workgroups
            )
        if workgroup_name == "":
            flash("Please choose a name for your Workgroup")
            return render_template(
                "create_workgroup.html", form=form, workgroups=workgroups
            )
        if not workgroup_name_delimiters_removed.isalnum():
            flash("Workgroup names cannot contain special characters!")
            return render_template(
                "create_workgroup.html", form=form, workgroups=workgroups
            )
        if info == "":
            flash("Please indicate why you need to set up a Workgroup")
            return render_template(
                "create_workgroup.html", form=form, workgroups=workgroups
            )

        # make sure workgroup of same name does not exist (add back in institution later)
        existing_workgroup_names = [
            x.name for x in db.session.query(models.WorkGroup).all()
        ]
        if [
            x
            for x in existing_workgroup_names
            if remove_spaces_and_dashes(x).lower()
            == workgroup_name_delimiters_removed.lower()
        ]:
            flash(
                "A Workgroup of this name already exists. Please choose a different name."
            )
            return render_template(
                "create_workgroup.html", form=form, workgroups=workgroups
            )

        existing_request = (
            db.session.query(models.WorkGroupRequest)
            .filter(func.lower(models.WorkGroupRequest.name) == workgroup_name.lower())
            .first()
        )
        if existing_request:
            flash(
                "This request has already been made. You will receive a notification when your request has been considered."
            )
            return render_template(
                "create_workgroup.html", form=form, workgroups=workgroups
            )

        # use the one default institution for now
        institution = (
            db.session.query(models.Institution)
            .filter(models.Institution.name == "Test User Institution")
            .first()
        )

        # create new workgroup and request for admin
        current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        new_workgroup = models.WorkGroup(
            name=workgroup_name,
            institution=institution.id,
            principal_investigator=[PI],
            approved=False,
            time_of_creation=current_time,
        )
        db.session.add(new_workgroup)
        db.session.commit()  # Have to commit to get the new_workgroup.id

        wg_request = models.WorkGroupRequest(
            name=workgroup_name,
            institution=institution.id,
            principal_investigator=PI.id,
            info=info,
            time=current_time,
            status="active",
            workgroup=new_workgroup.id,
        )
        db.session.add(wg_request)

        # send notification to admins
        admins = (
            db.session.query(models.User)
            .join(models.Role)
            .filter(models.Role.name == "Admin")
            .all()
        )
        for admin in admins:
            notification = models.Notification(
                person=admin.Person.id,
                type="New Workgroup Request",
                info="A PI has requested a new Workgroup. Go to the <a href=/admin_dashboard "
                'id="admin_dashboard">admin dashboard</a> to see this request.',
                wg_request="yes",
                time=datetime.now(),
                status="active",
                wg="",
                wb="",
            )
            db.session.add(notification)

            services.email.send_notification(admin.Person)
        db.session.commit()
        # flash success message and redirect to manage workgroup page
        flash(
            "Your Workgroup has been created. It will show as under moderation until it has been approved by the site admin."
            "You will receive a notification when your request has been approved."
        )
        return redirect(url_for("main.index"))
    return render_template("create_workgroup.html", form=form, workgroups=workgroups)
