from typing import List
from datetime import datetime

from flask_login import current_user
from sources import models, services
from sources.extensions import db


def send_and_commit(notification: models.Notification, person: models.Person) -> None:
    db.session.add(notification)
    db.session.commit()
    services.email.send_notification(person)


def deny_sm_status_request(person: models.Person) -> None:
    notification = models.Notification(
        person=person.id,
        type=f"Your Request to join {wg.name}",
        info=f"Your request to join Workgroup, {wg.name} has been denied.",
        time=datetime.now(),
        status="active",
    )

    send_and_commit(notification, person)


def deny_pi_status_request(person: models.Person) -> None:
    notification = models.Notification(
                person=person.id,
                type="Your Request to become a Principal Investigator",
                info=f"Your request to become a Principal Investigator in Workgroup, {wg.name}, has been denied.",
                time=datetime.now(),
                status="active",
            )

    send_and_commit(notification, person)


def deny_sr_status_request(person: models.Person) -> None:
    notification = models.Notification(
                person=person.id,
                type="Your Request to become a Senior Researcher",
                info=f"You request to become a Senior Researcher in Workgroup, {wg.name}, has been denied.",
                time=datetime.now(),
                status="active",
            )

    send_and_commit(notification, person)

def add_user_by_email_request(person: models.Person, workgroup: models.WorkGroup, role: str) -> models.Notification:
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
