import json
from typing import Optional

from flask import Response, flash, redirect, render_template, url_for
from flask_login import current_user, login_required
from sources import services
from sources.auxiliary import (
    get_notification_number,
    get_workgroups,
    security_member_workgroup,
)

from . import workgroup_bp


@workgroup_bp.route("/workgroup/<workgroup_selected>", methods=["GET", "POST"])
@workgroup_bp.route(
    "/workgroup/<workgroup_selected>/<workbook_selected>", methods=["GET", "POST"]
)
@login_required
@workgroup_bp.doc(security="sessionAuth")
def workgroup(
    workgroup_selected: str, workbook_selected: Optional[str] = None
) -> Response:
    """
    Workgroup page, shows all the workbooks in a workgroup and allows the user to select a workbook to view and
    see the reactions in that workbook. Must be a member of the workgroup.

    Args:
        workgroup_selected: workgroup to be selected
        workbook_selected: workbook to be selected

    Returns:
        flask.Response: renders the workgroup page
    """
    if not security_member_workgroup(workgroup_selected):
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))

    workgroups = get_workgroups()

    # finds workgroup object (needs institution later)
    workgroup_selected_obj = services.workgroup.from_name(workgroup_selected)
    approval_status = workgroup_selected_obj.approved
    notification_number = get_notification_number()

    user_type = services.workgroup.get_user_type(workgroup_selected, current_user)

    # get lists of workbooks, and the next reaction id for each workbooks
    workbooks = services.workbook.get_workbooks_from_user_group_combination(
        workgroup_selected
    )

    workbook_new_reaction_ids_dic = {}
    workbook_new_set_ids_dic = {}
    for workbook in workbooks:
        next_reaction_id = services.reaction.get_next_reaction_id_for_workbook(
            workbook.id
        )
        next_set_id = services.reaction_set.next_id_in_workbook(workbook.id)

        workbook_new_reaction_ids_dic[workbook.name] = next_reaction_id
        workbook_new_set_ids_dic[workbook.name] = next_set_id
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
        "work_structures/workgroup.html",
        workgroup_selected=workgroup_selected,
        workbooks=workbooks,
        workgroups=workgroups,
        user_type=user_type,
        notification_number=notification_number,
        workbook_selected=workbook_selected,
        approval_status=approval_status,
        workbook_next_reaction_ids=json.dumps(workbook_new_reaction_ids_dic),
        workbook_next_set_ids=json.dumps(workbook_new_set_ids_dic),
        new_reaction_id=new_reaction_id,
    )
