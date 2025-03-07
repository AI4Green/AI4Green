import base64
from datetime import date, datetime
from io import BytesIO

import qrcode
from dateutil.relativedelta import relativedelta
from flask import (
    Response,
    current_app,
    flash,
    jsonify,
    redirect,
    render_template,
    request,
    url_for,
)
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from PIL import Image
from sources import models, services
from sources.auxiliary import (
    get_all_workgroup_members,
    get_notification_number,
    get_workgroups,
    make_objects_inactive,
    security_member_workgroup,
    workgroup_dict,
)
from sources.decorators import (
    principal_investigator_required,
    workgroup_member_required,
)
from sources.extensions import db

from . import manage_workgroup_bp


@manage_workgroup_bp.route("/manage_workgroup/<workgroup>", methods=["GET", "POST"])
@manage_workgroup_bp.route(
    "/manage_workgroup/<workgroup>/<has_request>", methods=["GET", "POST"]
)
@login_required
@principal_investigator_required
@manage_workgroup_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator role in the specified workgroup.",
)
def manage_workgroup(workgroup: str, has_request: str = "no") -> Response:
    """
    Renders the page for managing a workgroup. This page allows the PI to view all members of the workgroup, remove
    members, change roles, and view requests to join the workgroup.

    Args:
        workgroup:
        has_request:

    Returns:
        flask.Response The rendered manage workgroup page for the specific workgroup
    """
    workgroups = get_workgroups()
    current_workgroup = workgroup
    notification_number = get_notification_number()
    # get all the pis, srs, and sms in the group
    wg = services.workgroup.from_name(current_workgroup)
    pi, sr, sm = get_all_workgroup_members(wg)

    # get requests for this workgroup
    user = (
        db.session.query(models.Person)
        .join(models.WorkGroup.principal_investigator)
        .filter(
            models.WorkGroup.principal_investigator.any(
                models.Person.user.has(email=current_user.email)
            )
        )
        .first()
    )

    requests = services.requests.get_active_in_workgroup_for_pi(user, wg)

    return render_template(
        "manage_workgroup.html",
        workgroups=workgroups,
        current_workgroup=current_workgroup,
        pi=pi,
        sr=sr,
        sm=sm,
        requests=requests,
        has_request=has_request,
        notification_number=notification_number,
        qr_code_expiration=date.today() + relativedelta(years=1),
    )


@manage_workgroup_bp.route(
    "/manage_workgroup/make_change/<workgroup>/<email>/<mode>/<current_status>",
    methods=["GET", "POST"],
)
@login_required
@principal_investigator_required
@manage_workgroup_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator role in the specified workgroup.",
)
def make_change_to_workgroup(
    workgroup: str, email: str, mode: str, current_status: str
) -> Response:
    """
    Enables a principal investigator to changes the role of a user in a workgroup or remove them from the workgroup.

    Args:
        workgroup: the workgroup name
        email: the email of the user to be changed or removed
        mode: the action to be taken, either remove or change role
        current_status: the current role of the user in the workgroup

    Returns:
        JSON message of success or failure
    """

    # find the person and workgroup objects that the change relates to
    wg = services.workgroup.from_name(workgroup)
    person = services.person.from_email(email)
    # mode is either remove user from workgroup or change role in workgroup
    if mode == "remove":
        if current_status == "pi" and len(wg.principal_investigator) < 2:
            return jsonify(
                {"feedback": "Cannot remove PI. Workgroups must have at least one PI."}
            )

        # uses the current_status variable as a dict key to determine the attribute being removed.
        getattr(person, workgroup_dict[current_status]["person_to_wg_attr"]).remove(wg)
        db.session.commit()
        return jsonify(
            {
                "feedback": f"{person.user.fullname} "
                f"has been removed from {wg.name} and associated Workbooks"
            }
        )
    else:
        old_role, new_role = mode.split("-to-")
        if old_role == "pi" and len(wg.principal_investigator) < 2:
            return jsonify(
                {"feedback": "Cannot remove PI. Workgroups must have at least one PI."}
            )
        # remove as old role then add as new one
        getattr(person, workgroup_dict[old_role]["person_to_wg_attr"]).remove(wg)
        getattr(person, workgroup_dict[new_role]["person_to_wg_attr"]).add(wg)
        db.session.commit()
        return jsonify(
            {
                "feedback": f"{person.user.fullname} has changed role from "
                f"{workgroup_dict[old_role]['display_string']} to "
                f"{workgroup_dict[new_role]['display_string']}"
            }
        )


# PI/SR status request page
@manage_workgroup_bp.route(
    "/PI_status_request/<user_type>/<new_role>/<workgroup>", methods=["GET", "POST"]
)
@manage_workgroup_bp.route(
    "/SR_status_request/<user_type>/<new_role>/<workgroup>", methods=["GET", "POST"]
)
@login_required
@workgroup_member_required
@manage_workgroup_bp.doc(
    security="sessionAuth",
    description="Requires login and membership in the specified workgroup.",
)
def status_request(user_type: str, new_role: str, workgroup: str) -> Response:
    """
    Allows a user to request a change in role within a workgroup. The request is sent to the workgroup

    Args:
        user_type: the current role of the user
        new_role: the requested new role of the user
        workgroup: the workgroup name

    Returns:
        Redirect to the home page with a message indicating the request has been made
    """
    # find all the PIs for the workgroup
    pis = (
        db.session.query(models.Person)
        .join(models.t_Person_WorkGroup)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
        .all()
    )
    wg = services.workgroup.from_name(workgroup)

    # find the person
    person = services.person.from_current_user_email()

    if services.notifications.duplicate_notification_check(
        [person], "New Workgroup Role Reassignment Request", "active", wg.name
    ):
        flash(
            "You have already submitted a request for this workgroup. You will receive a notification "
            "when your request has been considered."
        )
        return redirect(url_for("main.index"))

    # set up notification and request for each one
    for pi in pis:
        notification = models.Notification(
            person=pi.id,
            type="New Workgroup Role Reassignment Request",
            info=f"You have a new request from a member to change their role in Workgroup, {workgroup}, of which you are Principal Investigator.",
            time=datetime.now(),
            status="active",
            wg=wg.name,
        )
        db.session.add(notification)
        db.session.commit()  # needed to get the notification id.
        services.email.send_notification(pi)
        wg_request = models.WGStatusRequest(
            principal_investigator=pi.id,
            person=person.id,
            wg=wg.id,
            current_role=user_type,
            new_role=new_role,
            time=datetime.now(),
            status="active",
            notification=notification.id,
        )
        db.session.add(wg_request)
    db.session.commit()
    flash(
        "Your request has been made. You will receive a notification when your request has been considered."
    )
    return redirect(url_for("main.index"))


@manage_workgroup_bp.route(
    "/manage_workgroup/change_status_request/<workgroup>/<email>/<mode>/<decision>",
    methods=["GET", "POST"],
)
@login_required
@principal_investigator_required
@manage_workgroup_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator role in the specified workgroup.",
)
def change_status_from_request(
    workgroup: str, email: str, mode: str, decision: str
) -> Response:
    """
    Responds to a status request from a user to change their role in a workgroup. The PI can approve or deny the request
    and the user will be notified of the decision.

    Args:
        workgroup: the workgroup name
        email: the email of the user who requested the change
        mode: the mode of the request, either to-sm, sr-to-pi, or sm-to-pi
        decision: the decision of the PI, either approve or deny

    Returns:
        JSON message of success or failure
    """

    wg = services.workgroup.from_name(workgroup)

    # change request to inactive
    request_objs = services.requests.get_from_email_and_workgroup(email, wg)
    make_objects_inactive(request_objs)

    notification_objs = [x.Notification for x in request_objs]
    make_objects_inactive(notification_objs)

    person = services.person.from_email(email)

    if decision == "deny":
        # send new notification
        if mode == "to-sm":
            services.notifications.deny_sm_status_request(person, wg)
        elif mode == "sr-to-pi" or "sm-to-pi":
            services.notifications.deny_pi_status_request(person, wg)
        else:
            services.notifications.deny_sr_status_request(person, wg)

        return jsonify({"feedback": "This request has been denied!"})
    # decision is to approve user
    else:
        if len(mode.split("-to-")) > 1:
            old_role, new_role = mode.split("-to-")
            # remove old role
            getattr(person, workgroup_dict[old_role]["person_to_wg_attr"]).remove(wg)
            # set variables according to new role and workgroup
            new_role_display_string = workgroup_dict[new_role]["display_string"]
            request_type = f"Request to become a {new_role_display_string}"
            info = (
                f"Your request to become a {new_role_display_string} in Workgroup, {wg.name}, "
                f"has been approved!"
            )
            feedback = (
                f"{person.user.fullname} has changed role from {workgroup_dict[old_role]['display_string']} "
                f"to {workgroup_dict[new_role]['display_string']}"
            )
        else:
            # user is joining workgroup and has no old role
            new_role = "sm"
            request_type = "Request to join a Workgroup"
            info = f"Your request to join Workgroup {wg.name}, has been approved!"
            feedback = f"{person.user.fullname} has been added to {wg.name} as a Standard Member!"
        # add person to new role and send notification
        getattr(person, workgroup_dict[new_role]["person_to_wg_attr"]).add(wg)
        notification = models.Notification(
            person=person.id,
            type=request_type,
            info=info,
            time=datetime.now(),
            status="active",
        )
        services.notifications.send_and_commit(notification, person)
        return jsonify({"feedback": feedback})


@manage_workgroup_bp.route(
    "/manage_workgroup/go_to_workgroup/<workgroup>", methods=["GET", "POST"]
)
@login_required
@manage_workgroup_bp.doc(security="sessionAuth")
def go_to_workgroup(workgroup: str) -> Response:
    """
    Redirects the user to the workgroup management page from a notification
    Args:
        workgroup: the workgroup name

    Returns:
        Redirect to the workgroup management page
    """
    return redirect(
        url_for(
            "manage_workgroup.manage_workgroup", has_request="yes", workgroup=workgroup
        )
    )


@manage_workgroup_bp.route("/add_user_by_email", methods=["GET", "POST"])
@principal_investigator_required
@login_required
@manage_workgroup_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator role in the specified workgroup.",
)
def add_user_by_email(workgroup):
    """
    Allows a principal investigator to add a user to a workgroup by email. A request is sent to the user to join the
    workgroup.
    Args:
        workgroup: the workgroup name

    Returns:
        Redirect to the workgroup management page with a message indicating the request has been sent.
    """
    email = request.args.get("email")
    user_type = request.args.get("user_type")

    user = services.user.from_email(email)

    new_role_display_string = workgroup_dict[user_type]["display_string"]

    added_person = services.person.from_email(email)
    if not added_person:
        flash(f"User with email: {email} does not exist! Please try again.")
        return redirect(
            url_for(
                "manage_workgroup.manage_workgroup",
                workgroup=workgroup,
                has_request="no",
            )
        )

    wg = services.workgroup.from_name(workgroup)
    if not services.workgroup.get_user_type(workgroup, user):
        notification = services.notifications.add_user_by_email_request(
            added_person, wg, new_role_display_string
        )

        duplicates = services.requests.find_workgroup_duplicates_for_user(user, wg)

        if duplicates:
            flash(
                "This User's membership has already been requested for this Workgroup."
            )
            return redirect(
                url_for(
                    "manage_workgroup.manage_workgroup",
                    workgroup=workgroup,
                    has_request="no",
                )
            )

        else:
            services.notifications.send_and_commit(notification, added_person)

            wg_request = models.WGStatusRequest(
                principal_investigator=current_user.id,
                person=added_person.id,
                wg=wg.id,
                current_role="Non-Member",
                new_role=new_role_display_string,
                time=datetime.now(),
                status="active",
                notification=notification.id,
            )
            db.session.add(wg_request)
            db.session.commit()
            flash(
                f"A request has been send to User with email: {email} for this Workgroup!"
            )

    else:
        flash("This user is already a member of this Workgroup!")

    return redirect(
        url_for(
            "manage_workgroup.manage_workgroup", workgroup=workgroup, has_request="no"
        )
    )


@manage_workgroup_bp.route("/generate_qr_code/<workgroup>", methods=["GET", "POST"])
@principal_investigator_required
@login_required
@manage_workgroup_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator role in the specified workgroup.",
)
def generate_qr_code(workgroup=None):
    """
    Generates a QR code for a workgroup that allows a user to join the workgroup by scanning the code.

    Args:
        workgroup: the workgroup name

    Returns:
        flask.Response with a JSON message of the QR code image
    """
    token = services.email.get_encoded_token(31536000, {"workgroup": workgroup})
    url = f"https://{current_app.config['SERVER_NAME']}/qr_add_user/{token}"
    logo = Image.open("sources/static/img/favicon.ico")
    qr = qrcode.QRCode(
        version=1,
        box_size=10,
        border=4,
    )
    qr.add_data(url)
    qr.make(fit=True)

    img = qr.make_image(fill_color="#777778", back_color="white")
    pos = ((img.size[0] - logo.size[0]) // 2, (img.size[1] - logo.size[1]) // 2)
    img.paste(logo, pos)

    buffer = BytesIO()
    img.save(buffer, kind="PNG")
    qr_img = base64.b64encode(buffer.getvalue()).decode()

    return jsonify(qr_img)


@manage_workgroup_bp.route("/qr_add_user/<token>", methods=["GET", "POST"])
def add_user_by_qr(token=None):
    workgroup = services.email.verify_qr_code_for_add_user_token(token)
    if workgroup:
        return redirect(url_for("join_workgroup.join_workgroup", workgroup=workgroup))

    else:
        flash("This QR code has expired.")
        return redirect(url_for("main.index"))


@manage_workgroup_bp.route(
    "/change_workgroup_name/<workgroup>/<new_name>", methods=["POST"]
)
@principal_investigator_required
@login_required
@manage_workgroup_bp.doc(
    security="sessionAuth",
    description="Requires login and Principal Investigator role in the specified workgroup.",
)
def change_name(workgroup: str, new_name: str):
    """
    Calls a verify workgroup name function in /manage_workgroup/routes.py, verifies the output is the same as the user
    input. If so,the new name is committed to the database. If error, the feedback is returned with a 400 error.

    Args:
        workgroup_name: str, the current workgroup name
        new_name: str, the new name input by the user
    Returns:
        JSON message of either success or failure

    """
    verify_name = services.workgroup.verify_wg_name(workgroup, new_name)

    if verify_name == new_name:
        wg_object = services.workgroup.from_name(workgroup)
        wg_object.name = new_name
        db.session.commit()
        all_members = services.workgroup.list_all_members(wg_object)

        for member in all_members:
            services.notifications.new(member, workgroup, new_name)

        return (
            jsonify(
                {"message": "Workgroup name successfully changed", "new_name": new_name}
            ),
            200,
        )

    else:
        return jsonify({"error": verify_name}), 400
