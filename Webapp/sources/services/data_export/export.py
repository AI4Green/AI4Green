from typing import List

from flask import request
from flask_login import current_user
from sources import models, services
from sources.extensions import db


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
                self.workgroup.name, self.workbook.name, sort_crit="time"
            )
        self.workbooks = (
            set()
        )  # used to check the user has permission if reactions from multiple workbooks

    def check_permissions(self) -> str:
        """Turn reaction list into a workbook list. User must either be workbook member or PI"""
        # workbook_set = {}
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
        institution = services.workgroup.get_institution()
        reaction_id_list = [x.id for x in self.reaction_list]

        models.DataExportRequest.create(
            requestor=self.requestor.id,
            required_approvers=self.workgroup.principal_investigator,
            status="pending",
            reactions=reaction_id_list,
            workbooks=list(self.workbooks),
            workgroup=self.workgroup.id,
            institution=institution.id,
        )
