import json
from typing import Optional

from flask import Response, flash, redirect, render_template, url_for
from flask_login import current_user, login_required

from sources import models
from sources.auxiliary import (get_notification_number, get_workgroups,
                               security_member_workgroup)
from sources.extensions import db

from . import workgroup_bp


@workgroup_bp.route("/workgroup/<workgroup_selected>", methods=["GET", "POST"])
@workgroup_bp.route(
    "/workgroup/<workgroup_selected>/<workbook_selected>", methods=["GET", "POST"]
)
@login_required
def workgroup(
    workgroup_selected: str, workbook_selected: Optional[str] = None
) -> Response:
    # must be logged in and a member of the workgroup
    if not security_member_workgroup(workgroup_selected):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))
    workgroups = get_workgroups()

    # finds workgroup object (needs institution later)
    workgroup_selected_obj = (
        db.session.query(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_selected)
        .first()
    )
    approval_status = workgroup_selected_obj.approved
    notification_number = get_notification_number()

    pi = (
        db.session.query(models.User.email)
        .join(models.Person)
        .join(models.t_Person_WorkGroup)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_selected)
        .all()
    )
    if current_user.email in [user.email for user in pi]:
        user_type = "principal_investigator"

    sr = (
        db.session.query(models.User)
        .join(models.Person)
        .join(models.t_Person_WorkGroup_2)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_selected)
        .all()
    )
    if current_user.email in [user.email for user in sr]:
        user_type = "senior_researcher"

    sm = (
        db.session.query(models.User)
        .join(models.Person)
        .join(models.t_Person_WorkGroup_3)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_selected)
        .all()
    )
    if current_user.email in [user.email for user in sm]:
        user_type = "standard_member"

    # get lists of workbooks, and the next reaction id for each workbooks
    workbooks = (
        db.session.query(models.WorkBook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_selected)
        .join(models.t_Person_WorkBook)
        .join(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .all()
    )

    workbook_names = [wb.name for wb in workbooks]
    workbook_new_reaction_ids_dic = {}
    for workbook in workbooks:
        next_reaction_id = find_next_reaction_id(workbook)
        workbook_new_reaction_ids_dic[workbook.name] = next_reaction_id
    # select workbook with newest reaction in as active workbook to be selected by default
    if workbooks:
        workbook_object_list = (
            db.session.query(models.WorkBook)
            .filter(models.WorkBook.name.in_(workbook_names))
            .join(models.WorkGroup)
            .filter(models.WorkGroup.name == workgroup_selected)
            .all()
        )

        newest_reaction = (
            db.session.query(models.Reaction)
            .join(models.WorkBook)
            .filter(models.Reaction.workbooks.in_([x.id for x in workbook_object_list]))
            .order_by(models.Reaction.time_of_creation.desc())
            .first()
        )

        # if there is a newest reaction in the workbook collection, set to default and set the new reaction id
        if newest_reaction:
            workbook_selected_obj = newest_reaction.workbook
            # assign variable for workbook name and 3 letter abbreviation
            workbook_selected = workbook_selected_obj.name
        # if no reactions in any of the workbooks then take the workbook with the oldest primary key and make new reaction id
        else:
            workbook_selected = workbook_object_list[0].name
        # find the new reaction id of the active workbook
        new_reaction_id = workbook_new_reaction_ids_dic[workbook_selected]
    else:
        # if no workbooks assign variables here
        workbooks = "no_workbooks"
        new_reaction_id = None
    return render_template(
        "workgroup.html",
        workgroup_selected=workgroup_selected,
        workbooks=workbooks,
        workgroups=workgroups,
        user_type=user_type,
        notification_number=notification_number,
        workbook_selected=workbook_selected,
        approval_status=approval_status,
        workbook_next_reaction_ids=json.dumps(workbook_new_reaction_ids_dic),
        new_reaction_id=new_reaction_id,
    )


def find_next_reaction_id(workbook: models.WorkBook) -> str:
    workbook_abbreviation = workbook.abbreviation
    # find the newest reaction and then +1 to the id and return
    newest_reaction = (
        db.session.query(models.Reaction)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .order_by(models.Reaction.reaction_id.desc())
        .first()
    )
    if not newest_reaction:
        # if no reactions in workbook yet, then start with 001
        return f"{workbook_abbreviation}-001"
    most_recent_reaction_id = newest_reaction.reaction_id
    # take the number on the rhs of the reaction id, remove the 0s, convert to int, add 1, convert to str, add 0s
    new_reaction_id_number = str(
        int(most_recent_reaction_id.split("-")[-1].lstrip("0")) + 1
    ).zfill(3)
    return f"{workbook_abbreviation}-{new_reaction_id_number}"
