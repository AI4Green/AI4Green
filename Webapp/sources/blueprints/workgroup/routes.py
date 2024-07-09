import json
from typing import Optional

from flask import Response, flash, redirect, render_template, url_for
from flask_login import current_user, login_required

from sources import models, services
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
    workgroup_selected_obj = services.workgroup.from_name(workgroup_selected)
    approval_status = workgroup_selected_obj.approved
    notification_number = get_notification_number()

    user_type = services.workgroup.get_user_type(workgroup_selected)

    # get lists of workbooks, and the next reaction id for each workbooks
    workbooks = services.workbook.get_workbooks_from_user_group_combination(workgroup_selected)

    workbook_new_reaction_ids_dic = {}
    for workbook in workbooks:
        next_reaction_id = services.workbook.get_next_reaction_id_in_workbook(workbook)
        workbook_new_reaction_ids_dic[workbook.name] = next_reaction_id
    # select workbook with newest reaction in as active workbook to be selected by default
    if workbooks:

        newest_reaction = services.workbook.get_newest_reaction_in_workbooks(workbooks)

        # if there is a newest reaction in the workbook collection, set to default and set the new reaction id
        if newest_reaction:
            workbook_selected_obj = newest_reaction.workbook
            # assign variable for workbook name and 3 letter abbreviation
            workbook_selected = workbook_selected_obj.name

        # if no reactions in any of the workbooks then take the workbook with the oldest primary key and make new reaction id
        else:
            workbook_selected = workbooks[0].name

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
