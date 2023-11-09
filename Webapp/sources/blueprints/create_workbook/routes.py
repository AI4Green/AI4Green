from datetime import datetime

import pytz
from flask import redirect  # renders html templates
from flask import Response, flash, render_template, request, url_for
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from flask_wtf import FlaskForm
from sources import auxiliary, models
from sources.auxiliary import (
    get_notification_number,
    get_workgroups,
    remove_spaces_and_dashes,
    security_pi_sr_workgroup,
)
from sources.extensions import db
from sqlalchemy import func
from wtforms import StringField, SubmitField
from wtforms.validators import Length

from . import create_workbook_bp  # imports the blueprint of the dummy route


class CreateWorkbookForm(FlaskForm):
    """This form is used for creating a workbook"""

    workbook = StringField("Workbook", validators=[Length(min=4, max=140)])
    abbreviation = StringField(
        "3 letter abbreviation", validators=[Length(min=3, max=3)]
    )
    submit = SubmitField("Create Workbook", render_kw={"class": "btn btn-primary"})


@create_workbook_bp.route("/create_workbook/<workgroup>", methods=["GET", "POST"])
@login_required
def create_workbook(workgroup: str) -> Response:
    # must be logged in and a SR or PI of the workgroup
    if not security_pi_sr_workgroup(workgroup):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    """Creates a workbook using a FlaskForm. The book is created by a PI or SR and they must provide the name
    and workgroup the workbook should belong to"""
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    form = CreateWorkbookForm()

    wg = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .first()
    )
    # on valid submit
    if request.method == "POST" and form.validate_on_submit():
        # validates against special characters in workbook name
        workbook_name = auxiliary.sanitise_user_input(form.workbook.data)
        workbook_name_delimiters_removed = remove_spaces_and_dashes(workbook_name)
        if not workbook_name.replace(" ", "").replace("-", "").isalnum():
            flash("Workgroup names cannot contain special characters!")
            return render_template(
                "create_workgroup.html", form=form, workgroups=workgroups
            )

        # make sure workbook of same name does exist within workbook
        existing_workbook_names = (
            db.session.query(models.WorkBook)
            .join(models.WorkGroup)
            .filter(models.WorkGroup.id == wg.id)
            .all()
        )
        if [
            x
            for x in existing_workbook_names
            if remove_spaces_and_dashes(x.name).lower()
            == workbook_name_delimiters_removed.lower()
        ]:
            flash("A workbook of this name already exists!")
            return render_template(
                "create_workbook.html",
                form=form,
                workgroups=workgroups,
                notification_number=notification_number,
                workgroup=workgroup,
            )

        # validate workbook abbreviation
        workbook_abbreviation = auxiliary.sanitise_user_input(form.abbreviation.data)
        if not workbook_abbreviation.isalnum():
            flash("Workbook abbreviation must only consist of alphanumeric characters!")
            return render_template(
                "create_workbook.html",
                form=form,
                workgroups=workgroups,
                notification_number=notification_number,
                workgroup=workgroup,
            )
        existing_abbreviations = (
            db.session.query(models.WorkBook)
            .filter(
                func.lower(models.WorkBook.abbreviation)
                == workbook_abbreviation.lower()
            )
            .join(models.WorkGroup)
            .filter(models.WorkGroup.id == wg.id)
            .first()
        )
        if existing_abbreviations:
            flash("A workbook with this abbreviation already exists!")
            return render_template(
                "create_workbook.html",
                form=form,
                workgroups=workgroups,
                notification_number=notification_number,
                workgroup=workgroup,
            )

        # find user object
        user_obj = (
            db.session.query(models.Person)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .one()
        )

        # Create new workbook
        current_time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        workbook = models.WorkBook(
            name=workbook_name,
            abbreviation=workbook_abbreviation,
            group=wg.id,
            users=[user_obj],
            time_of_creation=current_time,
        )
        db.session.add(workbook)
        db.session.commit()

        wb_obj = (
            db.session.query(models.WorkBook)
            .filter(models.WorkBook.name == workbook_name)
            .join(models.WorkGroup)
            .filter(models.WorkGroup.name == workgroup)
            .first()
        )
        # flash success message and redirect to manage workgroup page
        flash("Workbook has been created!")
        return redirect(
            url_for(
                "manage_workbook.manage_workbook",
                workgroup=workgroup,
                has_request="no",
                workbook=wb_obj.id,
            )
        )
    else:
        return render_template(
            "create_workbook.html",
            form=form,
            workgroups=workgroups,
            notification_number=notification_number,
            workgroup=workgroup,
        )
