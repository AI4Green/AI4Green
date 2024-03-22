import threading
import time
from datetime import datetime
from typing import List, Literal

import pytz
from flask import flash, request
from flask_login import current_user
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


class RequestStatus:
    def __init__(self, data_export_request: models.DataExportRequest):
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
        Updates the status of a data export request to 'REJECTED' and notifies the requestor of rejection.
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
            f"{self.data_export_request.data_format} was rejected.",
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
            print("APPROVED!")
            self.data_export_request.status = "APPROVED"
            db.session.commit()


class RequestLinkVerification:
    """Verifies access to data export request response page."""

    def __init__(self, token):
        self.token = token
        self.approver = None
        self.data_export_request = None
        self.status = ""

    def verify_request_link(self):
        """Checks user is a required approver for this request and that it has not already been approved."""
        (
            self.approver,
            self.data_export_request,
        ) = services.email.verify_data_export_request_token(self.token)
        if not (
            self.approver == current_user
            and services.person.from_current_user_email()
            in self.data_export_request.required_approvers
        ):
            flash("Data export request expired")
            return "failed"
        # redirect if user has already approved.
        if self._check_if_already_responded() is True:
            flash("You have already responded to this data export request")
            return "failed"
        return "success"

    def _check_if_already_responded(self) -> bool:
        """Returns True if the user has already responded to this data export request"""
        return (
            db.session.query(models.data_export_request_approvers.c.responded)
            .filter(
                models.data_export_request_approvers.c.data_export_request_id
                == self.data_export_request.id
            )
            .filter(
                models.data_export_request_approvers.c.person_id == self.approver.id
            )
            .first()
        )[0]


class InitiateDataExport:
    def __init__(self, data_export_request: models.DataExportRequest):
        self.data_export_request = data_export_request

    def initiate(self):
        export_function = self.get_export_function()
        export_function()
        print(export_function.__name__)
        print("work")
        time.sleep(20)
        print("still working")

    def get_export_function(self):
        export_function_dict = {
            "RDF": self.make_rdf_export,
            "CSV": self.make_csv_export,
            "JSON": self.make_json_export,
            "SURF": self.make_surf_export,
            "PDF": self.make_pdf_export,
            "ELN": self.make_eln_export,
        }
        return export_function_dict[self.data_export_request.data_format.value]

    def make_rdf_export(self):
        # container_name = hash(self.data_export_request)
        # for reaction in self.data_export_request.reactions:
        #     rdf = services.data_export.reaction_data_file.ReactionDataFile(
        #         reaction,
        #         reaction.reaction_id,
        #     )
        pass

    def make_csv_export(self):
        pass

    def make_json_export(self):
        pass

    def make_surf_export(self):
        pass

    def make_pdf_export(self):
        pass

    def make_eln_export(self):
        pass


class MakeZip:
    """Makes a zip file of the reaction data files"""

    def __init__(self, workbook, workgroup):
        self.workbook = services.workbook.get_workbook_from_group_book_name_combination(
            workbook, workgroup
        )
        self.workgroup = workgroup
        self.requestor = current_user
        # self.make_zip()

    #
    #     def make_zip(self):
    #         """Temporarily saves a zip file to the blob service"""
    #
    #         reactions = services.reaction.list_active_in_workbook(
    #             self.workbook.name, self.workgroup
    #         )
    #         rdf_list = []
    #
    #         for reaction in reactions:
    #             ReactionDataFile(
    #                 reaction,
    #             )
