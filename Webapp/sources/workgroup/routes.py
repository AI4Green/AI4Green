from flask import render_template, redirect, url_for, flash
from sources.workgroup import workgroup_bp
from flask_login import login_required
from pony.orm import select, desc
from sources import db
from flask_login import current_user
from sources.auxiliary import get_workgroups, get_notification_number, security_member_workgroup, sanitise_user_input
import json



@workgroup_bp.route('/workgroup/<workgroup_selected>', methods=['GET', 'POST'])
@workgroup_bp.route('/workgroup/<workgroup_selected>/<workbook_selected>', methods=['GET', 'POST'])
@login_required
def workgroup(workgroup_selected, workbook_selected=None):
    # must be logged in and a member of the workgroup
    if not security_member_workgroup(workgroup_selected):
        flash("You do not have permission to view this page")
        return redirect(url_for('main.index'))
    workgroups = get_workgroups()
    # finds workgroup object (needs institution later)
    workgroup_selected_obj = select(x for x in db.WorkGroup if x.name == workgroup_selected).first()
    approval_status = workgroup_selected_obj.approved
    notification_number = get_notification_number()
    pi = select(x.principal_investigator.user.email for x in db.WorkGroup if x.name == workgroup_selected)[:]
    if current_user.email in pi:
        user_type = "principal_investigator"
    sr = select(x.senior_researcher.user.email for x in db.WorkGroup if x.name == workgroup_selected)[:]
    if current_user.email in sr:
        user_type = "senior_researcher"
    sm = select(x.standard_member.user.email for x in db.WorkGroup if x.name == workgroup_selected)[:]
    if current_user.email in sm:
        user_type = "standard_member"
    # get lists of workbooks, and the next reaction id for each workbooks
    workbooks = []
    for workbook in select(u for u in db.WorkBook if workgroup_selected == u.group.name and current_user.email in u.users.user.email):
        workbooks.append(workbook)
    workbook_names = [wb.name for wb in workbooks]
    workbook_new_reaction_ids_dic = {}
    for workbook in workbooks:
        next_reaction_id = find_next_reaction_id(workbook)
        workbook_new_reaction_ids_dic[workbook.name] = next_reaction_id
    # select workbook with newest reaction in as active workbook to be selected by default
    if workbooks:
        workbook_object_list = select(wb for wb in db.WorkBook if wb.group.name == workgroup_selected and wb.name in workbook_names).order_by(
            lambda x: x.id)[:]
        newest_reaction = select(rxn for rxn in db.Reaction if rxn.workbooks in workbook_object_list).order_by(
            lambda r: desc(r.time_of_creation)).first()
        # if there is a newest reaction in the workbook collection, set to default and set the new reaction id
        if newest_reaction:
            workbook_selected_obj = newest_reaction.workbooks
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
    return render_template("workgroup.html", workgroup_selected=workgroup_selected, workbooks=workbooks,
                           workgroups=workgroups, user_type=user_type, notification_number=notification_number,
                           workbook_selected=workbook_selected, approval_status=approval_status,
                           workbook_next_reaction_ids=json.dumps(workbook_new_reaction_ids_dic), new_reaction_id=new_reaction_id)


def find_next_reaction_id(workbook):
    workbook_abbreviation = workbook.abbreviation
    # find the newest reaction and then +1 to the id and return
    newest_reaction = select(rxn for rxn in db.Reaction if rxn.workbooks == workbook).order_by(
        lambda r: desc(r.reaction_id)).first()
    if not newest_reaction:
        # if no reactions in workbook yet, then start with 001
        return workbook_abbreviation + '-001'
    most_recent_reaction_id = newest_reaction.reaction_id
    # take the number on the rhs of the reaction id, remove the 0s, convert to int, add 1, convert to str, add 0s
    new_reaction_id_number = str(int(most_recent_reaction_id.split('-')[-1].lstrip('0')) + 1).zfill(3)
    new_reaction_id = workbook_abbreviation + '-' + new_reaction_id_number
    return new_reaction_id
