from datetime import datetime
from typing import List, Optional

from flask_login import current_user
from sources import models, services
from sources.extensions import db


def send_and_commit(notification: models.Notification, person: models.Person) -> None:
    """
    Sends notification email to user and commits notification to database
    Args:
        notification: models.Notification, notification to commit
        person: models.Person, who to send the email to

    Returns:
        None
    """
    db.session.add(notification)
    db.session.commit()
    services.email.send_notification(person)


def deny_sm_status_request(person: models.Person, workgroup: models.WorkGroup) -> None:
    """
    Creates notification denying a users request for joining a workgroup
    Args:
        person: models.Person, Person to send the notification to
        workgroup: models.WorkGroup, Workgroup they are trying to join

    Returns:
        None
    """
    notification = models.Notification(
        person=person.id,
        type=f"Your Request to join {workgroup.name}",
        info=f"Your request to join Workgroup, {workgroup.name} has been denied.",
        time=datetime.now(),
        status="active",
    )

    send_and_commit(notification, person)


def deny_pi_status_request(person: models.Person, workgroup: models.WorkGroup) -> None:
    """
    Creates notification denying a users request for PI status in workgroup
    Args:
        person: models.Person, Person to send the notification to
        workgroup: models.WorkGroup, Workgroup they are trying to join

    Returns:
        None
    """
    notification = models.Notification(
        person=person.id,
        type="Your Request to become a Principal Investigator",
        info=f"Your request to become a Principal Investigator in Workgroup, {workgroup.name}, has been denied.",
        time=datetime.now(),
        status="active",
    )

    send_and_commit(notification, person)


def deny_sr_status_request(person: models.Person, workgroup: models.WorkGroup) -> None:
    """
    Creates notification denying a users request for senior researcher status
    Args:
        person: models.Person, Person to send the notification to
        workgroup: models.WorkGroup, Workgroup they are trying to join

    Returns:
        None
    """
    notification = models.Notification(
        person=person.id,
        type="Your Request to become a Senior Researcher",
        info=f"You request to become a Senior Researcher in Workgroup, {workgroup.name}, has been denied.",
        time=datetime.now(),
        status="active",
    )

    send_and_commit(notification, person)


def add_user_by_email_request(
    person: models.Person, workgroup: models.WorkGroup, role: str
) -> models.Notification:
    """
    Sends request to user when they have been added to a workgroup by email.
    Args:
        person: models.Person, who to send the notification to
        workgroup: models.WorkGroup, the workgroup they will be added to
        role: str, role in the workgroup once added. either principal_investigator, senior_researcher or standard_member

    Returns:
        models.Notification, the notification object to be added to the database.
        Note that this is not committed or sent in this function.
    """
    return models.Notification(
        person=person.id,
        type="You Have Been Added to a Workgroup",
        info="A Principal Investigator has added you to the Workgroup, "
        + workgroup.name
        + ", as a "
        + role
        + ". Please respond below.",
        time=datetime.now(),
        status="active",
        wg=workgroup.name,
        wb="",
        wg_request="added",
    )


def duplicate_notification_check(
    people_ls: List[models.Person],
    request_type: str,
    status,
    WG: str,
    WB: models.WorkBook = None,
) -> bool:
    """
    Check if the request has already been made - if someone already has received your request they will be
    removed from the list of people to be sent the request
    """
    refined_people_ls = []
    for person in people_ls:
        if WB:
            check = (
                db.session.query(models.Notification)
                .filter(models.Notification.wb == WB and models.Notification.wg == WG)
                .filter(models.Notification.status == status)
                .join(models.WBStatusRequest)
                .join(models.Person, models.WBStatusRequest.person == models.Person.id)
                .filter(models.Person.id == person.id)
                .first()
            )
        else:
            if request_type == "New Workgroup Membership Request":
                # check for person applying to join workgroup - 1 PI may have many requests from different users
                check = (
                    db.session.query(models.Notification)
                    .filter(models.Notification.wg == WG)
                    .join(models.WGStatusRequest)
                    .join(
                        models.Person, models.WGStatusRequest.person == models.Person.id
                    )
                    .filter(models.Person == person.id)
                    .first()
                )
            elif request_type == "New Workgroup Role Reassignment Request":
                check = (
                    db.session.query(models.Notification)
                    .filter(models.Notification.wg == WG)
                    .filter(models.Notification.status == status)
                    .join(models.WGStatusRequest)
                    .join(
                        models.Person, models.WGStatusRequest.person == models.Person.id
                    )
                    .filter(models.Person.id == person.id)
                    .first()
                )
            else:
                check = (
                    db.session.query(models.Notification)
                    .filter(models.Notification.wg == WG)
                    .filter(models.Notification.status == status)
                    .filter(models.Notification.type == request_type)
                    .join(models.WGStatusRequest)
                    .join(
                        models.Person, models.WGStatusRequest.person == models.Person.id
                    )
                    .filter(models.Person.id == person.id)
                    .first()
                )

        if check is None:
            refined_people_ls.append(person)
    return not refined_people_ls


def new(person: models.Person, old_name: str, new_name: str) -> None:
    """
    Sends notification on AI4Green to user and commits notification to database
    Args:
        person: models.Person, person to whose sent information regarding workgroup name change
        old_name: str, the old name of the workgroup
        new_name: str, the new name of the workgroup
    Returns:
        None
    """
    notification = models.Notification(
        person=person.id,
        type="Your Workgroup has been renamed",
        info=f"The Workgroup name: {old_name} has been changed to {new_name}",
        time=datetime.now(),
        status="active",
    )

    db.session.add(notification)
    db.session.commit()
