from datetime import datetime
from typing import List, Literal

import pytz
from flask import request
from sources import models, services
from sources.extensions import db
from sqlalchemy import update


class NewRequest:
    def __init__(self):
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
            the permission result is either accepted or denied
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

    def _user_in_workbook_or_pi(self, workbook) -> bool:
        return (
            self.requestor in workbook.users
            or self.requestor in workbook.workgroup.principal_investigator
        )

    def _get_reaction_list(self) -> List[models.Reaction]:
        """Either get customised reaction list from request or all reactions in a workbook."""
        if request.json["reactionIDList"]:
            reaction_list = [
                services.reaction.get_from_reaction_id_and_workbook_id(
                    rxn_id, self.workbook.id
                )
                for rxn_id in request.json["reactionIDList"]
            ]
        else:
            reaction_list = services.reaction.list_active_in_workbook(
                self.workbook.name, self.workgroup.name, sort_crit="time"
            )
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
                type="Data Export Request",  # todo add link to info
                info=f"The following user: {self.requestor.user.email} has requested to export data from workbooks: "
                f"{', '.join(workbook_names)} in {self.data_format} format.<br>"
                f"To respond please follow the link sent to your email account. This request will expire after 7 days.",
                time=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
                status="active",
            )
            services.email.send_data_export_approval_request(
                principal_investigator.user, self.data_export_request
            )


def request_rejected(data_export_request: models.DataExportRequest):
    """
    Updates the status of a data export request to 'REJECTED' and notifies the requestor of rejection.

    Args:
        data_export_request (models.DataExportRequest): The data export request object.
    """
    data_export_request.status = "REJECTED"
    db.session.commit()
    # notify requestor of rejection
    workbooks = ", ".join([wb.name for wb in data_export_request.workbooks])
    models.Notification.create(
        person=data_export_request.requestor,
        type="Data Export Request Rejected",
        info=f"Your request to export {len(data_export_request.reactions)} reactions from {workbooks} in "
        f"{data_export_request.data_format} was rejected.",
        time=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
        status="active",
    )


def request_accepted(data_export_request: models.DataExportRequest):
    """
    Accepts a data export request, updates the approval status of the current user's approval
    in the data export request approvers association table, and updates the status of the data export request.

    Args:
        data_export_request - The data export request object for which the current user's approval is being updated.
    """
    approver = services.person.from_current_user_email()
    # use update query because it is an association table. Normal query returns immutable Row.
    update_query = (
        update(models.data_export_request_approvers)
        .where(
            models.data_export_request_approvers.c.data_export_request_id
            == data_export_request.id
        )
        .where(models.data_export_request_approvers.c.person_id == approver.id)
        .values(approved=True)
    )
    db.session.execute(update_query)
    db.session.commit()


def update_request_status(data_export_request: models.DataExportRequest):
    """
    Updates the status of a data export request based on the approval decisions of the approvers.

    Args:
       data_export_request (models.DataExportRequest): The data export request object.
    """
    # TODO update this to work with association table
    data_export_approvers = models.data_export_request_approvers.query.filter_by(
        data_export_request_id=data_export_request.id
    ).first()
    approval_decisions = []
    for approver in data_export_approvers:
        approval_decisions.append(
            models.data_export_request_approvers.query.filter_by(
                data_export_request_id=data_export_request.id, person_id=approver.id
            )
            .first()
            .approved
        )

    if all(approval_decisions is True):
        data_export_request.status = "APPROVED"
        db.session.commit()
        # initiate data export release with threads then notify requestor after successful data export generation.
