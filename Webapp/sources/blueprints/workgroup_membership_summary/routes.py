from flask import Response, render_template  # renders html templates
from flask_login import \
    login_required  # protects a view function against anonymous users
from flask_wtf import FlaskForm
from wtforms import SelectField, SubmitField

from sources.auxiliary import (get_notification_number, get_workbooks,
                               get_workgroups)

from . import \
    workgroup_membership_summary_bp  # imports the blueprint of the dummy route


class SelectWorkgroupForm(FlaskForm):
    workgroups = SelectField("Workgroup")
    submit = SubmitField("Find Workbooks")


# put routes here
@workgroup_membership_summary_bp.route(
    "/workgroup_membership_summary", methods=["GET", "POST"]
)
@login_required
def workgroup_membership_summary() -> Response:
    # must be logged in
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    # create the form
    form = SelectWorkgroupForm()
    workgroups.insert(0, "-select-")
    form.workgroups.choices = workgroups
    # on valid submit
    if form.validate_on_submit() and form.workgroups.data != "-select-":
        workbooks = get_workbooks(form.workgroups.data)
        if not workbooks:
            workbooks.append("You are not a member of any workbooks in this workgroup")
        return render_template(
            "workgroup_membership_summary.html",
            form=form,
            workgroup=form.workgroups.data,
            workbooks=workbooks,
            workgroups=workgroups,
            notification_number=notification_number,
        )
    return render_template(
        "workgroup_membership_summary.html",
        form=form,
        workgroups=workgroups,
        notification_number=notification_number,
    )
