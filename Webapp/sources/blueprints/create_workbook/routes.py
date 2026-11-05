from datetime import datetime

import pytz
from flask import Response, flash, redirect, render_template, request, url_for
from flask_login import current_user, login_required
from flask_wtf import FlaskForm
from sources import auxiliary, models, services
from sources.auxiliary import get_notification_number, get_workgroups
from sources.decorators import principal_investigator_or_senior_researcher_required
from sources.extensions import db
from sqlalchemy import func
from wtforms import StringField, SubmitField
from wtforms.validators import Length

from . import create_workbook_bp


class CreateWorkbookForm(FlaskForm):
    """This form is used for creating a workbook"""

    workbook = StringField("Workbook", validators=[Length(min=4, max=140)])
    abbreviation = StringField(
        "3 letter abbreviation", validators=[Length(min=3, max=3)]
    )
    submit = SubmitField("Create Workbook", render_kw={"class": "btn btn-primary"})


@create_workbook_bp.route("/create_workbook/<workgroup>", methods=["GET", "POST"])
@login_required
@principal_investigator_or_senior_researcher_required
@create_workbook_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator or Senior Researcher role in the specified workgroup.",
)
def create_workbook(workgroup: str) -> Response:
    """Creates a workbook using a FlaskForm. The book is created by a principal investigator
     or senior researcher and they must provide the name and workgroup the workbook should belong to

    Args:
        workgroup (str): The name of the workgroup the workbook should belong to

    Returns:
        flask.Response: The rendered template for creating a workbook
        or a redirect to the manage workbook page for the new workbook
    """
    workgroup_name = workgroup  # make variable names consistent with other files
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    form = CreateWorkbookForm()

    workgroup = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )
    # on valid submit
    if request.method == "POST" and form.validate_on_submit():
        # validates against special characters in workbook name
        workbook_name = auxiliary.sanitise_user_input(form.workbook.data)
        workbook_name_delimiters_removed = services.utils.remove_spaces_and_dashes(
            workbook_name
        )
        if not workbook_name.replace(" ", "").replace("-", "").isalnum():
            flash("Workgroup names cannot contain special characters!")
            return render_template(
                "work_structures/create_workgroup.html",
                form=form,
                workgroups=workgroups,
            )

        # make sure workbook of same name does exist within workbook
        existing_workbook_names = (
            db.session.query(models.WorkBook)
            .join(models.WorkGroup)
            .filter(models.WorkGroup.id == workgroup.id)
            .all()
        )
        if [
            x
            for x in existing_workbook_names
            if services.utils.remove_spaces_and_dashes(x.name).lower()
            == workbook_name_delimiters_removed.lower()
        ]:
            flash("A workbook of this name already exists!")
            return render_template(
                "work_structures/create_workbook.html",
                form=form,
                workgroups=workgroups,
                notification_number=notification_number,
                workgroup=workgroup_name,
            )

        # validate workbook abbreviation
        workbook_abbreviation = auxiliary.sanitise_user_input(form.abbreviation.data)
        if not workbook_abbreviation.isalnum():
            flash("Workbook abbreviation must only consist of alphanumeric characters!")
            return render_template(
                "work_structures/create_workbook.html",
                form=form,
                workgroups=workgroups,
                notification_number=notification_number,
                workgroup=workgroup_name,
            )
        existing_abbreviations = (
            db.session.query(models.WorkBook)
            .filter(
                func.lower(models.WorkBook.abbreviation)
                == workbook_abbreviation.lower()
            )
            .join(models.WorkGroup)
            .filter(models.WorkGroup.id == workgroup.id)
            .first()
        )
        if existing_abbreviations:
            flash("A workbook with this abbreviation already exists!")
            return render_template(
                "work_structures/create_workbook.html",
                form=form,
                workgroups=workgroups,
                notification_number=notification_number,
                workgroup=workgroup_name,
            )

        # find person object
        person = (
            db.session.query(models.Person)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .one()
        )

        # Create new workbook
        current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        new_workbook = models.WorkBook(
            name=workbook_name,
            abbreviation=workbook_abbreviation,
            group=workgroup.id,
            users=[person],
            time_of_creation=current_time,
        )
        db.session.add(new_workbook)
        db.session.commit()

        workbook = (
            db.session.query(models.WorkBook)
            .filter(models.WorkBook.name == workbook_name)
            .join(models.WorkGroup)
            .filter(models.WorkGroup.name == workgroup_name)
            .first()
        )

        # record access change
        message = services.data_access_history.DataAccessMessage(
            person.user.fullname,
            person.user.email,
            workgroup.id,
            "No Access",
            "Access",
            datetime.now().strftime("%Y-%m-%d"),
            workbook.id,
        )
        services.data_access_history.send_message(message)

        # flash success message and redirect to manage workgroup_name page
        flash("Workbook has been created!")
        return redirect(
            url_for(
                "manage_workbook.manage_workbook",
                workgroup=workgroup_name,
                has_request="no",
                workbook=workbook.id,
            )
        )
    else:
        return render_template(
            "work_structures/create_workbook.html",
            form=form,
            workgroups=workgroups,
            notification_number=notification_number,
            workgroup=workgroup_name,
        )
