from flask import render_template, request, flash, redirect, url_for  # renders html templates
from flask_login import login_required  # protects a view function against anonymous users
from sources.delete_profile import delete_profile_bp  # imports the blueprint of the route
from sources.auxiliary import get_workgroups, get_notification_number
from flask_login import current_user
from sources import db
from pony.orm import delete, select
from datetime import datetime
from sources.email_methods import send_notification_email


@delete_profile_bp.route('/confirm_delete_profile', methods=['GET', 'POST'])
@login_required
def confirm_delete_profile():
    # must be logged in
    if request.method == "POST":
        # remove from all workgroups and workbooks
        wgs_pi = select(x for x in db.WorkGroup if current_user.email in x.principal_investigator.user.email)[:]
        wgs_sr = select(x for x in db.WorkGroup if current_user.email in x.senior_researcher.user.email)[:]
        wgs_sm = select(x for x in db.WorkGroup if current_user.email in x.standard_member.user.email)[:]
        for wg in wgs_pi:
            person = select(x for x in wg.principal_investigator if x.user.email == current_user.email).first()
            db.Person[person.id].workgroup_principal_investigator.remove(db.WorkGroup[wg.id])
            workbooks_to_remove = list((select(b.id for b in db.WorkBook if b.group.id == wg.id)[:]))
            for workbook in workbooks_to_remove:
                db.Person[person.id].workbook_user.remove(db.WorkBook[workbook])
            # if wg orphaned send notification to admins
            if len(db.WorkGroup[wg.id].principal_investigator) == 0:
                admins = select(u for u in db.User if u.role.name == "Admin")[:]
                for admin in admins:
                    admin_person = select(u for u in db.Person if u.user.email == admin.email).first()
                    db.Notification(person=admin_person, type="A Workgroup has been archived",
                                    info="A PI has deleted their profile and Workgroup, " + wg.name +
                                         ", has been archived.", time=datetime.now(), status="active")
                    send_notification_email(admin_person)
                # archive workgroup by removing all members
                for p in select(u for u in db.WorkGroup[wg.id].principal_investigator)[:]:
                    db.Person[p.id].workgroup_principal_investigator.remove(db.WorkGroup[wg.id])
                for p in select(u for u in db.WorkGroup[wg.id].senior_researcher)[:]:
                    db.Person[p.id].workgroup_senior_researcher.remove(db.WorkGroup[wg.id])
                for p in select(u for u in db.WorkGroup[wg.id].standard_member)[:]:
                    db.Person[p.id].workgroup_standard_member.remove(db.WorkGroup[wg.id])
        for wg in wgs_sr:
            person = select(x for x in wg.senior_researcher if x.user.email == current_user.email).first()
            db.Person[person.id].workgroup_senior_researcher.remove(db.WorkGroup[wg.id])
            workbooks_to_remove = list((select(b.id for b in db.WorkBook if b.group.id == wg.id)[:]))
            for workbook in workbooks_to_remove:
                db.Person[person.id].workbook_user.remove(db.WorkBook[workbook])
        for wg in wgs_sm:
            person = select(x for x in wg.standard_member if x.user.email == current_user.email).first()
            db.Person[person.id].workgroup_standard_member.remove(db.WorkGroup[wg.id])
            workbooks_to_remove = list((select(b.id for b in db.WorkBook if b.group.id == wg.id)[:]))
            for workbook in workbooks_to_remove:
                db.Person[person.id].workbook_user.remove(db.WorkBook[workbook])
        # delete any requests and associated notifications
        delete(x.notification for x in db.WGStatusRequest if x.person.user.email == current_user.email)
        delete(x.notification for x in db.WBStatusRequest if x.person.user.email == current_user.email)
        delete(x for x in db.WGStatusRequest if x.person.user.email == current_user.email)
        delete(x for x in db.WBStatusRequest if x.person.user.email == current_user.email)
        delete(x for x in db.WorkGroup_request if x.principal_investigator.user.email == current_user.email)
        # delete profile
        delete(p for p in db.User if p == current_user)
        flash("Your profile has been deleted!")
        return redirect(url_for("auth.logout"))
    wgs_pi = select(x for x in db.WorkGroup if current_user.email in x.principal_investigator.user.email and
                    len(x.principal_investigator) == 1)[:]
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template("delete_confirm.html", workgroups=workgroups, notification_number=notification_number,
                           wgs_pi=wgs_pi)


@delete_profile_bp.route('/delete_profile/reassign/<wg_name>', methods=['GET', 'POST'])
@login_required
def delete_profile_reassign(wg_name):
    return redirect(url_for("manage_workgroup.manage_workgroup", workgroup=wg_name))
