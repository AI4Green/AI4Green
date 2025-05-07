from datetime import datetime, timedelta
from typing import Dict, List, Union

import pytz
from flask import request
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook
from sources.extensions import db
from sqlalchemy import func, update


def get_from_name_and_workbook_id(name: str, workbook_id: int) -> models.Reaction:
    """
    Retrieves a reaction based on its name and workbook ID.

    Args:
        name (str): The name of the reaction.
        workbook_id (int): The ID of the workbook to which the reaction belongs.

    Returns:
        models.Reaction: The reaction object matching the name and workbook ID.
    """
    return (
        db.session.query(models.Reaction)
        .filter(func.lower(models.Reaction.name) == name.lower())
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook_id)
        .first()
    )


def get_current_from_request() -> models.Reaction:
    """
    Gets the current reaction for a request from the frontend. Either using request.form or request.json
    Returns:
        The reaction which matches the details of the request
    """
    if request.form:
        reaction = get_current_from_request_form()
    elif request.json:
        reaction = get_current_from_request_json()
    return reaction


def get_current_from_request_form() -> models.Reaction:
    """
    Gets the current reaction using the request.form variable
    Returns:
        Reaction that corresponds to data in request.form
    """
    reaction_id = str(request.form["reactionID"])
    workgroup_name = str(request.form["workgroup"])
    workbook_name = str(request.form["workbook"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    return (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )


def get_current_from_request_json() -> models.Reaction:
    """
    Gets the current reaction using request.json
    Returns:
        Reaction that corresponds to data in request.json
    """
    reaction_id = str(request.json["reactionID"])
    workgroup_name = str(request.json["workgroup"])
    workbook_name = str(request.json["workbook"])
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    return (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )


def get_current_from_request_args() -> models.Reaction:
    """
    Gets the current reaction using the request.args variable
    Returns:
        Reaction that corresponds to data in request.args
    """
    reaction_id = str(request.args.get("reactionID"))
    workgroup_name = str(request.args.get("workgroup"))
    workbook_name = str(request.args.get("workbook"))
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)
    return (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook.id)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup_name)
        .first()
    )


def list_recent() -> List[models.Reaction]:
    """
    Gets a list of reactions created in the past 28 days. For the admin_dashboard

    Returns:
         List of all reactions from past 28 days
    """
    cut_off_date = datetime.now(pytz.timezone("Europe/London")).replace(
        tzinfo=None
    ) - timedelta(days=28)
    return (
        (
            db.session.query(models.Reaction).filter(
                models.Reaction.time_of_creation > cut_off_date
            )
        )
        .order_by(models.Reaction.time_of_creation.desc())
        .all()
    )


def count() -> int:
    """
    Gets the number of reactions in the database

    Returns:
        The number of reactions in the database
    """
    return db.session.query(models.Reaction).count()


def get_from_reaction_id_and_workbook_id(
    reaction_id: str, workbook_id: int
) -> models.Reaction:
    """
    Gets the reaction from the reaction_id and workbook id

    Args:
        reaction_id - in format WB1-001
        workbook - The workbook the reaction belongs to

    Returns:
        models.Reaction: The reaction object matching the reaction ID and workbook ID.
    """
    return (
        db.session.query(models.Reaction)
        .filter(models.Reaction.reaction_id == reaction_id)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook_id)
        .first()
    )


def list_active_in_workbook(
    workbook: str, workgroup: str, sort_crit: str = "AZ"
) -> List[models.Reaction]:
    """
    Gets the active reactions for a workbook. Active means the reaction has not been deleted/archived.

    Args:
        workbook (str): The name of the workbook.
        workgroup (str): The name of the workgroup.
        sort_crit (str, optional): The sorting criteria for the reaction list.
            Defaults to 'AZ' for alphabetical sorting.

    Returns:
        List[models.Reaction]: A list of active reactions in the specified workbook and workgroup,
        sorted based on the specified criteria ('AZ' for alphabetical, 'time' for time of creation).
    """

    query = (
        db.session.query(models.Reaction)
        .filter(models.Reaction.status == "active")
        .join(models.WorkBook)
        .filter(models.WorkBook.name == workbook)
        .join(models.WorkGroup)
        .filter(models.WorkGroup.name == workgroup)
    )
    if sort_crit == "time":
        reaction_list = query.order_by(models.Reaction.time_of_creation.desc()).all()
    elif sort_crit == "AZ":
        reaction_list = query.order_by(models.Reaction.name.asc()).all()
    return reaction_list


def make_reaction_image_list(reaction_list: List[models.Reaction]) -> List[str]:
    """
    gets a list of reaction images from db
    Args:
        reaction_list: list of Reactions objects that we are making scheme images for

    Returns:
        A list of reaction images as strings
    """
    reaction_image_list = []
    # get images from db
    for reaction in reaction_list:
        reaction_image = reaction.reaction_image
        # if there is no image use a blank string to maintain correct list length
        if reaction_image:
            reaction_image_list.append(reaction_image)
        else:
            reaction_image_list.append("")
    return reaction_image_list


def to_dict(reaction_list: List[models.Reaction]) -> List[Dict]:
    """
    Converts a reaction list to a dictionary used to render template: '_saved_reactions.html'

    Args:
        reaction_list - list of reactions as objects
        sort_crit - criteria to sort reactions.

    Returns:
        A List of dictionaries with the reaction data required to render the _saved_reactions.html template

    """

    reactions = []
    for idx, reaction in enumerate(reaction_list):
        # for each reaction get the relevant info and shorten description if it's long
        description = reaction.description
        if reaction.creator_person.user:
            creator_email = reaction.creator_person.user.email
            creator_username = reaction.creator_person.user.username
        else:
            creator_email = "unknown"
            creator_username = "a deleted profile"
        if len(description) > 250:
            description = description[0:249] + "..."
        reaction_details = {
            "html_id": idx + 1,
            "name": reaction.name,
            "description": description,
            "time_of_creation": str(reaction.time_of_creation),
            "time_of_update": str(reaction.time_of_update),
            "reaction_smiles": reaction.reaction_smiles,
            "reaction_table_data": reaction.reaction_table_data,
            "summary_table_data": reaction.summary_table_data,
            "workgroup": reaction.WorkBook.WorkGroup.name,
            "workbook": reaction.WorkBook.name,
            "completion_status": reaction.complete,
            "reaction_id": reaction.reaction_id,
            "creator_email": creator_email,
            "creator_username": creator_username,
            "addenda": reaction.addenda,
        }
        reactions.append(reaction_details)
    return reactions


def add_addendum(
    reaction: models.Reaction, reaction_note_text: str, author: models.Person
) -> models.ReactionNote:
    """
    Add a new addendum to a reaction.

    Args:
        reaction - The reaction object to which the addendum is added.
        reaction_note_text - The text content of the addendum.
        author - The Person object for the author of the addendum

    Returns:
        The newly created reactionNote

    """
    new_addendum = models.ReactionNote(
        text=reaction_note_text,
        time_of_creation=datetime.now(),
        author=author.id,
        reaction=reaction.id,
    )
    db.session.add(new_addendum)
    db.session.commit()
    return new_addendum


def get_addenda(reaction: models.Reaction) -> List[models.ReactionNote]:
    """
    Get all addenda for a reaction
    Args:
        reaction: The reaction object for which addenda are retrieved

    Returns:
        List of ReactionNote (addenda) objects

    """
    return (
        db.session.query(models.ReactionNote)
        .join(models.Reaction)
        .filter(models.Reaction.id == reaction.id)
        .all()
    )


def most_recent_in_workbook(workbook_id: int) -> models.Reaction:
    """
    Retrieves the most recent reaction in a given workbook.

    Args:
        workbook_id (int): The ID of the workbook.

    Returns:
        models.Reaction: The most recent reaction object in the workbook.
    """
    return (
        db.session.query(models.Reaction)
        .join(models.WorkBook)
        .filter(models.WorkBook.id == workbook_id)
        .order_by(models.Reaction.reaction_id.desc())
        .first()
    )


def get_next_reaction_id_for_workbook(workbook_id: int) -> str:
    """
    Generates the next reaction ID for a given workbook in format WB1-001

    Args:
        workbook_id (int): The ID of the workbook.

    Returns:
        str: The next reaction ID for the workbook.
    """
    workbook_obj = services.workbook.get(workbook_id)
    workbook_abbreviation = workbook_obj.abbreviation
    # find the newest reaction and then +1 to the id and return
    newest_reaction = most_recent_in_workbook(workbook_id)
    if not newest_reaction:
        # if no reactions in workbook yet, then start with 001
        return workbook_abbreviation + "-001"
    most_recent_reaction_id = newest_reaction.reaction_id
    # take the number on the rhs of the reaction id, remove the 0s, convert to int, add 1, convert to str, add 0s
    new_reaction_id_number = str(
        int(most_recent_reaction_id.split("-")[-1].lstrip("0")) + 1
    ).zfill(3)
    new_reaction_id = workbook_abbreviation + "-" + new_reaction_id_number
    return new_reaction_id


def add(
    name: str,
    reaction_id: str,
    creator: models.Person,
    workbook_id: int,
    reaction_table: Dict[str, any],
    summary_table: Dict[str, any],
    reaction_smiles: str = "",
) -> models.Reaction:
    """
    Adds a reaction to the database.

    Args:
        name (str): The name of the reaction.
        reaction_id (str): The generated id for the reaction in format WB1-001
        creator (models.Person): The creator of the reaction.
        workbook_id (int): The ID of the workbook to which the reaction belongs.
        reaction_table (Dict[str, any]): Data for the reaction table.
        summary_table (Dict[str, any]): Data for the summary table.
        reaction_smiles Optional(str): The SMILES representation of the reaction.

    Returns:
        models.Reaction: The newly added reaction object.
    """
    reaction = models.Reaction(
        name=name,
        reaction_id=reaction_id,
        creator=creator.id,
        workbooks=workbook_id,
        status="active",
        complete="not complete",
        reaction_smiles=reaction_smiles,
        reaction_table_data=reaction_table,
        summary_table_data=summary_table,
    )
    db.session.add(reaction)
    db.session.commit()
    return reaction


class NewReactionApprovalRequest:
    """
    Class to handle creation of new reaction approval requests
    """

    def __init__(self):
        """Creates an instance of NewReactionApprovalRequest from request.json"""
        self.requestor = services.user.person_from_current_user()
        self.workgroup = services.workgroup.from_name(request.json.get("workgroup"))
        self.workbook = services.workbook.get_workbook_from_group_book_name_combination(
            self.workgroup.name, request.json.get("workbook")
        )
        self.reaction = get_current_from_request_json()
        self.reaction_approval_request = None

    def submit_request(self):
        """Save request to database and notify and email approvers."""
        self._save_to_database()
        self._notify_approvers()

    def _save_to_database(self):
        """Save the data export request to the database"""
        institution = services.workgroup.get_institution()

        self.reaction_approval_request = models.ReactionApprovalRequest.create(
            requestor=self.requestor.id,
            required_approvers=self.workgroup.principal_investigator,
            status="PENDING",
            reaction=self.reaction.id,
            workgroup=self.workgroup.id,
            workbook=self.workbook.id,
            institution=institution.id,
        )

    def _notify_approvers(self):
        """
        Send a notification and an email to each principal investigator of the reaction workgroup
        """
        for principal_investigator in self.workgroup.principal_investigator:
            models.Notification.create(
                person=principal_investigator.id,
                type="Reaction Approval Request",
                info=self._message_content(),
                time=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
                status="active",
            )
            services.email.send_reaction_approval_request(
                principal_investigator.user,
                self.reaction_approval_request,
                self.workgroup,
                self.workbook,
                self.reaction,
            )

    def _message_content(self):
        return f"""
                <p>A user has requested your review for their reaction!</p>
                <table style='border-collapse: collapse; width: 100%; max-width: 600px; font-family: Arial, sans-serif;'>
                  <tr style='background-color: #f2f2f2;'>
                    <th style='text-align: center; padding: 8px; border: 1px solid #ddd;'>User</th>
                    <th style='text-align: center; padding: 8px; border: 1px solid #ddd;'>Workgroup</th>
                    <th style='text-align: center; padding: 8px; border: 1px solid #ddd;'>Workbook</th>
                    <th style='text-align: center; padding: 8px; border: 1px solid #ddd;'>Reaction Title</th>
                    <th style='text-align: center; padding: 8px; border: 1px solid #ddd;'>Reaction ID</th>
                  </tr>
                  <tr>
                    <td style='padding: 8px; border: 1px solid #ddd;'>{self.requestor.user.email}</td>
                    <td style='padding: 8px; border: 1px solid #ddd;'>{self.workgroup.name}</td>
                    <td style='padding: 8px; border: 1px solid #ddd;'>{self.workbook.name}</td>
                    <td style='padding: 8px; border: 1px solid #ddd;'>{self.reaction.name}</td>
                    <td style='padding: 8px; border: 1px solid #ddd;'>{self.reaction.reaction_id}</td>
                  </tr>
                </table>
                <br>
                <button class='btn btn-primary'>Review Reaction</button>
                """


class ReactionApprovalRequestStatus:
    """Class to update the request status when an approver approves or rejects a reaction"""

    def __init__(self, reaction_approval_request: models.ReactionApprovalRequest):
        """Creates an instance of the RequestStatus class."""
        self.reaction_approval_request = reaction_approval_request
        self.approver = services.person.from_current_user_email()

    def _update_query(self, approved: bool):
        """
        Docstring for update_query.
        """
        # use update query because it is an association table. Normal query returns immutable Row.
        update_query = (
            update(models.reaction_approval_request_approvers)
            .where(
                models.reaction_approval_request_approvers.c.reaction_approval_request_id
                == self.reaction_approval_request.id
            )
            .where(
                models.reaction_approval_request_approvers.c.person_id
                == self.approver.id
            )
            .values(approved=approved, responded=True)
        )
        db.session.execute(update_query)

    def _update_request_status(self, status: str, comments: Union[str, None] = None):
        """
        status is ENUM, how to document?
        """

        self.reaction_approval_request.status = status
        self.reaction_approval_request.reviewed_by = self.approver.id
        self.reaction_approval_request.time_of_review = datetime.now()
        if comments is not None:
            self.reaction_approval_request.comments = comments

        db.session.commit()

    def approve(self):
        """
        Approves a reaction approval request, updates the approval status of the current user's approval
        in the data export request approvers association table, and updates the status of the reaction approval request.
        """
        self._update_query(True)
        self._update_request_status("APPROVED")

    def reject(self, comments):
        """
        Updates the status of a reaction approval request to as responded to (bool) and 'REJECTED' and notifies the requestor.
        """
        self._update_query(False)
        self._update_request_status("REJECTED")

        # notify requestor of rejection
        models.Notification.create(
            person=self.reaction_approval_request.requestor,
            type="Reaction Rejected",
            info=f"Your reaction {self.reaction_approval_request.reaction.reaction_id}  in workbook "
            f"{self.reaction_approval_request.workbook.name} has been rejected.",
            time=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
            status="active",
        )

    def suggest_changes(self, comments):
        """
        Updates the status of a reaction approval request to as responded to (bool) and 'CHANGES_REQUESTED' and notifies the requestor.
        """
        self._update_query(False)
        self._update_request_status("CHANGES_REQUESTED", comments=comments)
