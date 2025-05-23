import html

from flask import Response, jsonify, render_template, request
from flask_login import current_user, login_required
from sources import models
from sources.auxiliary import get_notification_number, get_workgroups
from sources.extensions import db

from . import notifications_bp


@notifications_bp.route("/notifications", methods=["GET", "POST"])
@login_required
@notifications_bp.doc(security="sessionAuth")
def notifications() -> Response:
    """
    Renders the notification page to display the notifications for the current user

    Returns:
        flask.Response: The notifications page with all the notifications for the current user
    """
    workgroups = get_workgroups()
    notifications_obj = (
        db.session.query(models.Notification)
        .filter(models.Notification.status != "inactive")
        .join(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .order_by(models.Notification.time.desc())
        .all()
    )
    for notification in notifications_obj:
        notification.status = "read"
    db.session.commit()
    notification_number = get_notification_number()
    return render_template(
        "account_management/notifications.html",
        notifications=notifications_obj,
        workgroups=workgroups,
        notification_number=notification_number,
    )


@notifications_bp.route("/archive_notification", methods=["GET", "POST"])
@login_required
@notifications_bp.doc(security="sessionAuth")
def archive_notification() -> Response:
    """
    Archive a notification

    Returns:
        flask.Response: A JSON response indicating whether the notification was successfully archived

    """
    type_str = str(request.form["type"])
    info_str = str(request.form["info"])
    time = str(request.form["time"])
    info_str = html.unescape(info_str)
    notifications_obj = (
        db.session.query(models.Notification)
        .filter(models.Notification.type == type_str)
        .filter(models.Notification.info == info_str)
        .filter(models.Notification.time == time)
        .join(models.Person)
        .join(models.User)
        .filter(models.User.email == current_user.email)
        .order_by(models.Notification.time.desc())
        .first()
    )
    notifications_obj.status = "inactive"
    db.session.commit()
    return jsonify({"feedback": "Success"})
