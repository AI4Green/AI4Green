from flask import render_template, request, jsonify
from flask_login import login_required, current_user
from sources.notifications import notifications_bp
from pony.orm import select, desc
from sources import db
from sources.auxiliary import get_workgroups, get_notification_number
from html.parser import HTMLParser


@notifications_bp.route('/notifications', methods=['GET', 'POST'])
@login_required
def notifications():
    # must be logged in
    workgroups = get_workgroups()
    notifications_obj = select(x for x in db.Notification if x.person.user.email == current_user.email and x.status !=
                               "inactive").order_by(desc(db.Notification.time))[:]
    for notification in notifications_obj:
        notification.status = "read"
    notification_number = get_notification_number()
    return render_template("notifications.html", notifications=notifications_obj, workgroups=workgroups,
                           notification_number=notification_number)


@notifications_bp.route('/archive_notification', methods=['GET', 'POST'])
@login_required
def archive_notification():
    # must be logged in
    type_str = str(request.form['type'])
    info_str = str(request.form['info'])
    time = str(request.form['time'])
    parser = HTMLParser()
    info_str = parser.unescape(info_str)
    notifications_obj = select(x for x in db.Notification if x.person.user.email == current_user.email and x.type ==
                               type_str and x.info == info_str and str(x.time) == time).first()
    notifications_obj.status = "inactive"
    return jsonify({'feedback': 'Success'})
