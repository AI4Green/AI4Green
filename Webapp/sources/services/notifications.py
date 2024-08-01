from typing import List
from datetime import datetime

from flask_login import current_user
from sources import models, services
from sources.extensions import db


def send_and_commit(notification: models.Notification, person: models.Person) -> None:
    """
    sends notification email to user and commits notification to database
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

def add_user_by_email_request(person: models.Person, workgroup: models.WorkGroup, role: str) -> models.Notification:
    """
    sends request to user when they have been added to a workgroup by email.
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
