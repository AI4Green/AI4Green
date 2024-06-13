from datetime import datetime
from typing import List, Literal

import pytz
from flask import flash, request
from flask_login import current_user
from sources import db, models, services
from sqlalchemy import update


class NewRequest:
    """Class to handle the creation of new data export requests"""

    def __init__(self):
        """Creates an instance of the NewRequest class using from the request json."""
        self.requestor = services.user.person_from_current_user()
        self.data_format = request.json["exportFormat"]
        self.invalid_workbooks = None
        self.workgroup = services.workgroup.from_name(request.json["workgroup"])
        self.workbook = services.workbook.get_workbook_from_group_book_name_combination(
            self.workgroup.name, request.json["workbook"]
        )
        # if the user has customised the list of reactions to export use that else export all in workbook.
        self.reaction_list = self._get_reaction_list()
        # used to check the user has permission if reactions from multiple workbooks
        self.workbooks = set()
        self.data_export_request = None

    def check_permissions(
        self,
    ) -> Literal["permission denied", "permission accepted", "no data to export"]:
        """
        Turn reaction list into a workbook list. User must either be workbook member or PI
        Returns:
            the permission result is either accepted or denied or no data
        """
        # if there are no reactions we cannot export
        if len(self.reaction_list) == 0:
            return "no data to export"
        # get each workbook that we are exporting data from
        self._get_workbooks()
        # check if the user is either a PI of the workgroup or a workbook member for each workbook.
        for workbook in self.workbooks:
            if self._user_in_workbook_or_pi(workbook):
                continue
            else:
                return "permission denied"
        return "permission accepted"

    def submit_request(self):
        """Save request to database and notify and email approvers."""
        self._save_to_database()
        self._notify_approvers()

    def _get_workbooks(self):
        """Adds each workbook we are exporting data from to the self.workbooks set."""
        for reaction in self.reaction_list:
            self.workbooks.add(
                services.workbook.get_workbook_from_group_book_name_combination(
                    reaction.workbook.WorkGroup.name, reaction.workbook.name
                )
            )

    def _user_in_workbook_or_pi(self, workbook: models.workbook) -> bool:
        """Returns True if the user is either in the workbook or is the PI of the workgroup."""
        return (
            self.requestor in workbook.users
            or self.requestor in workbook.workgroup.principal_investigator
        )

    def _get_reaction_list(self) -> List[models.Reaction]:
        """Get customised reaction list from request"""
        reaction_list = [
            services.reaction.get_from_reaction_id_and_workbook_id(
                rxn_id, self.workbook.id
            )
            for rxn_id in request.json["reactionIDList"]
        ]
        return reaction_list

    def _save_to_database(self):
        """Save the data export request to the database"""
        institution = services.workgroup.get_institution()

        self.data_export_request = models.DataExportRequest.create(
            data_format=self.data_format,
            requestor=self.requestor.id,
            required_approvers=self.workgroup.principal_investigator,
            status="PENDING",
            reactions=self.reaction_list,
            workbooks=list(self.workbooks),
            workgroup=self.workgroup.id,
            institution=institution.id,
        )

    def _notify_approvers(self):
        """
        Send a notification and an email to each principal investigator of the workgroup the data is being exported from
        """
        workbook_names = [wb.name for wb in self.workbooks]
        for principal_investigator in self.workgroup.principal_investigator:
            models.Notification.create(
                person=principal_investigator.id,
                type="Data Export Request",
                info=f"The following user: {self.requestor.user.email} has requested to export data from workbooks: "
                f"{', '.join(workbook_names)} in {self.data_format} format.<br>"
                f"To respond please follow the link sent to your email account. This request will expire after 7 days.",
                time=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
                status="active",
            )
            services.email.send_data_export_approval_request(
                principal_investigator.user, self.data_export_request
            )


class RequestStatus:
    """Class to update the request status when an approver accepts or denies a request"""

    def __init__(self, data_export_request: models.DataExportRequest):
        """Creates an instance of the RequestStatus class."""
        self.data_export_request = data_export_request

    def accept(self):
        """
        Accepts a data export request, updates the approval status of the current user's approval
        in the data export request approvers association table, and updates the status of the data export request.
        """
        approver = services.person.from_current_user_email()
        # use update query because it is an association table. Normal query returns immutable Row.
        update_query = (
            update(models.data_export_request_approvers)
            .where(
                models.data_export_request_approvers.c.data_export_request_id
                == self.data_export_request.id
            )
            .where(models.data_export_request_approvers.c.person_id == approver.id)
            .values(approved=True, responded=True)
        )
        db.session.execute(update_query)
        db.session.commit()

    def deny(self):
        """
        Updates the status of a data export request to as responded to (bool) and 'REJECTED' and notifies the requestor.
        """
        approver = services.person.from_current_user_email()
        update_query = (
            update(models.data_export_request_approvers)
            .where(
                models.data_export_request_approvers.c.data_export_request_id
                == self.data_export_request.id
            )
            .where(models.data_export_request_approvers.c.person_id == approver.id)
            .values(approved=False, responded=True)
        )
        db.session.execute(update_query)
        db.session.commit()
        self.data_export_request.status = "REJECTED"
        db.session.commit()
        # notify requestor of rejection
        workbooks = ", ".join([wb.name for wb in self.data_export_request.workbooks])
        models.Notification.create(
            person=self.data_export_request.requestor,
            type="Data Export Request Rejected",
            info=f"Your request to export {len(self.data_export_request.reactions)} reactions from {workbooks} in "
            f"{self.data_export_request.data_format} was rejected by a principal investigator in your workgroup.",
            time=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
            status="active",
        )

    def update_status(self):
        """
        Updates the status attribute of a data export request if all approvers have approved.
        """
        # get all the required approvers and check if they have all approved.
        data_export_approvers = (
            db.session.query(models.data_export_request_approvers.c.approved)
            .filter(
                models.data_export_request_approvers.c.data_export_request_id
                == self.data_export_request.id
            )
            .all()
        )
        # get the results out of the row objects.
        approval_results = [x[0] for x in data_export_approvers]

        if all(x is True for x in approval_results):
            self.data_export_request.status = "APPROVED"
            db.session.commit()


class RequestLinkVerification:
    """Verifies access to data export request pages when using links containing tokens."""

    def __init__(self, token: str):
        self.token = token
        self.person = None
        self.data_export_request = None
        self.status = ""

    def verify_requestor_link(self) -> str:
        """We check the requestor link to access the download which is given once the request has been approved"""
        (
            self.person,
            self.data_export_request,
        ) = services.email.verify_data_export_token(self.token)
        # redirect if requestor does not match current user and request has not been approved.
        if not (
            self.person
            and self.data_export_request
            and self.person == current_user
            and services.person.from_current_user_email()
            == self.data_export_request.requestor_person
            and self.data_export_request.status.value == "APPROVED"
        ):
            flash("Data export request expired")
            return "failed"
        return "success"

    def verify_approver_link(self) -> str:
        """Checks user is a required approver for this request and that it has not already been approved."""
        (
            self.person,
            self.data_export_request,
        ) = services.email.verify_data_export_token(self.token)
        if not (
            self.person
            and self.data_export_request
            and self.person == current_user
            and services.person.from_current_user_email()
            in self.data_export_request.required_approvers
        ):
            flash("Data export request expired")
            return "failed"
        # redirect if user has already approved.
        if False:  # TODO debugging dev code
            # if self._check_if_already_responded() is True:
            flash("You have already responded to this data export request")
            return "failed"
        return "success"

    def _check_if_already_responded(self) -> bool:
        """Returns True if the user has already responded to this data export request"""
        approver_responded = (
            db.session.query(models.data_export_request_approvers)
            .filter(
                models.data_export_request_approvers.c.data_export_request_id
                == self.data_export_request.id
            )
            .filter(
                models.data_export_request_approvers.c.person_id
                == self.person.Person.id
            )
            .first()
        )
        return approver_responded[-1]
