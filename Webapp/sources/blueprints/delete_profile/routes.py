import contextlib
from datetime import datetime

from flask import redirect  # renders html templates
from flask import Response, flash, render_template, request, url_for
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from sources import models, services
from sources.auxiliary import get_notification_number, get_workgroups
from sources.extensions import db

from . import delete_profile_bp


@delete_profile_bp.route("/delete_profile", methods=["GET", "POST"])
@login_required
@delete_profile_bp.doc(security="sessionAuth")
def delete_profile() -> Response:
    """
    Renders the delete user profile page

    Returns:
        flask.Response: The rendered template for deleting a user profile
    """
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "account_management/delete_profile.html",
        workgroups=workgroups,
        notification_number=notification_number,
    )


@delete_profile_bp.route("/confirm_delete_profile", methods=["GET", "POST"])
@login_required
@delete_profile_bp.doc(security="sessionAuth")
def confirm_delete_profile() -> Response:
    """
    Confirms the deletion of a user profile

    Returns:
        flask.Response: The rendered template for confirming the deletion of a user.
    """
    if request.method == "POST":
        # remove from all workgroups and workbooks
        wgs_pi = (
            db.session.query(models.WorkGroup)
            .join(models.t_Person_WorkGroup)
            .join(models.Person)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .all()
        )
        wgs_sr = (
            db.session.query(models.WorkGroup)
            .join(models.t_Person_WorkGroup_2)
            .join(models.Person)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .all()
        )
        wgs_sm = (
            db.session.query(models.WorkGroup)
            .join(models.t_Person_WorkGroup_3)
            .join(models.Person)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .all()
        )

        person = (
            db.session.query(models.Person)
            .join(models.User)
            .filter(models.User.email == current_user.email)
            .first()
        )
        for wg in wgs_pi:
            person.workgroup_principal_investigator.remove(wg)

            workbooks_to_remove = (
                db.session.query(models.WorkBook)
                .join(models.WorkGroup)
                .filter(models.WorkGroup.id == wg.id)
                .all()
            )
            for workbook in workbooks_to_remove:
                person.workbook_user.remove(workbook)

            # if wg orphaned send notification to admins
            if len(wg.principal_investigator) == 0:
                admins = (
                    db.session.query(models.Person)
                    .join(models.User)
                    .join(models.Role)
                    .filter(models.Role.name == "Admin")
                    .all()
                )
                for admin in admins:
                    notification = models.Notification(
                        person=admin.id,
                        type="A Workgroup has been archived",
                        info=f"A PI has deleted their profile and Workgroup, {wg.name}, has been archived.",
                        time=datetime.now(),
                        status="active",
                    )
                    db.session.add(notification)
                    services.email_services.send_notification(admin)
                # archive workgroup by removing all members
                for p in wg.senior_researcher:
                    p.workgroup_senior_researcher.remove(wg)
                for p in wg.standard_member:
                    p.workgroup_senior_researcher.remove(wg)
        for wg in wgs_sr:
            person.workgroup_senior_researcher.remove(wg)
            for workbook in wg.WorkBook:
                with contextlib.suppress(ValueError):
                    person.workbook_user.remove(workbook)
        for wg in wgs_sm:
            person.workgroup_standard_member.remove(wg)
            for workbook in wg.WorkBook:
                with contextlib.suppress(ValueError):
                    person.workbook_user.remove(workbook)

        # delete profile and associated objects
        user = (
            db.session.query(models.User)
            .filter(models.User.email == current_user.email)
            .first()
        )

        objects_to_delete = [
            person.person_notification,
            person.person_status_change,
            person.person_status_change_pi,
            person.person_status_change_wb,
            person.person_status_change_wb_sr_pi,
        ]

        for obj_lists in objects_to_delete:
            for obj in obj_lists:
                db.session.delete(obj)

        db.session.delete(user)
        db.session.commit()

        # record access change
        for wg in wgs_pi + wgs_sr + wgs_sm:
            message = services.data_access_history.DataAccessMessage(
                person.id,
                wg.id,
                old_role="Access",
                new_role="No Access",
                date=datetime.now().strftime("%Y-%m-%d"),
            )
            services.data_access_history.send_message(message)

        flash("Your profile has been deleted!")
        return redirect(url_for("auth.logout"))

    wgs_pi = (
        db.session.query(models.WorkGroup)
        .join(models.t_Person_WorkGroup)
        .join(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .all()
    )
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    return render_template(
        "account_management/delete_confirm.html",
        workgroups=workgroups,
        notification_number=notification_number,
        wgs_pi=wgs_pi,
    )


@delete_profile_bp.route("/delete_profile/reassign/<wg_name>", methods=["GET", "POST"])
@login_required
@delete_profile_bp.doc(security="sessionAuth")
def delete_profile_reassign(wg_name: str) -> Response:
    """
    Redirects the user to the workgroup to reassign their principal investigator role before deleting their profile

    Args:
        wg_name (str): The name of the workgroup the user is reassigning their role in.

    Returns:
        flask.Response: A redirect to the workgroup management page to reassign their principal investigator role
    """
    return redirect(url_for("manage_workgroup.manage_workgroup", workgroup=wg_name))
