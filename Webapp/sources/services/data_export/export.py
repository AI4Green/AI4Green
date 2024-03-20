from datetime import datetime

import pytz
from flask import request
from sources import models, services


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
        if request.json["reactionIDList"]:
            self.reaction_list = [
                services.reaction.get_from_reaction_id_and_workbook_id(
                    rxn_id, self.workbook.id
                )
                for rxn_id in request.json["reactionIDList"]
            ]
        else:
            self.reaction_list = services.reaction.list_active_in_workbook(
                self.workbook.name, self.workgroup.name, sort_crit="time"
            )
        self.workbooks = (
            set()
        )  # used to check the user has permission if reactions from multiple workbooks

    def check_permissions(self) -> str:
        """Turn reaction list into a workbook list. User must either be workbook member or PI"""
        # workbook_set = {}
        if len(self.reaction_list) == 0:
            return "no data to export"
        for reaction in self.reaction_list:
            print(type(reaction))
            print(reaction.workbook.WorkGroup.name, reaction.workbook.name)
            self.workbooks.add(
                services.workbook.get_workbook_from_group_book_name_combination(
                    reaction.workbook.WorkGroup.name, reaction.workbook.name
                )
            )
        for workbook in self.workbooks:
            if (
                self.requestor in workbook.users
                or self.requestor in workbook.workgroup.principal_investigator
            ):
                continue
            else:
                return "permission denied"
        return "permission accepted"

    def submit_request(self):
        """Save request to database and notify and email approvers."""
        self._save_to_database()
        self._notify_approvers()

    def _save_to_database(self):
        institution = services.workgroup.get_institution()

        models.DataExportRequest.create(
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
        workbook_names = [wb.name for wb in self.workbooks]
        for principal_investigator in self.workgroup.principal_investigator:
            models.Notification.create(
                person=principal_investigator.id,
                type="Your request to export data",  # todo add link to info
                info=f"The following user: {self.requestor.user.fullname} has requested to export data from workbooks: "
                f"{', '.join(workbook_names)} in {self.data_format} format.\n"
                f"This request will expire after 7 days.",
                time=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
                status="active",
            )
            services.email.send_notification(principal_investigator)
