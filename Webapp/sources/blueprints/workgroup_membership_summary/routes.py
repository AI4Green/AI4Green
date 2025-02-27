from flask import Response, render_template
from flask_login import login_required
from flask_wtf import FlaskForm
from sources.auxiliary import get_notification_number, get_workbooks, get_workgroups
from wtforms import SelectField, SubmitField

from . import workgroup_membership_summary_bp


class SelectWorkgroupForm(FlaskForm):
    workgroups = SelectField("Workgroup")
    submit = SubmitField("Find Workbooks")


@workgroup_membership_summary_bp.route(
    "/workgroup_membership_summary", methods=["GET", "POST"]
)
@login_required
def workgroup_membership_summary() -> Response:
    """
    Workgroup membership summary page, all the workgroups the user belongs to and can find the workbooks
    in each workgroup they are a member of.

    Returns:
        flask.Response: renders the workgroup membership summary page
    """
    workgroups = get_workgroups() or ["-select-"]
    # create the form
    form = SelectWorkgroupForm()
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
        )
    return render_template(
        "workgroup_membership_summary.html",
        form=form,
        workgroups=workgroups,
    )
