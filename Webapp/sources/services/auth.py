from typing import Dict

from flask import abort, request
from flask_login import current_user
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook


def reaction_files(
    permission_level: str,
    request_source: str = "user",
    file_object_for_deletion: models.ReactionDataFile = None,
):
    """
    Authenticates user to either view or edit the reaction.
    Args:
        permission_level: Takes value of 'edit' or 'view_only'
        request_source: The origin of the request, either 'user' or 'server'
        file_object_for_deletion: If deleting a file object, we pass that object as an argument.
    """
    if permission_level == "view_only":
        view_files()
    elif permission_level == "edit":
        reaction = services.reaction.get_current_from_request()
        edit_reaction(
            reaction,
            request_source,
            file_attachment=True,
            file_object_for_deletion=file_object_for_deletion,
        )


def view_files():
    """Authenticates user as a workbook member or aborts. Gets the workgroup_name, workbook_name, and workbook."""
    if request.method == "GET":
        workgroup_name = request.args.get("workgroup")
        workbook_name = request.args.get("workbook")
    else:
        workgroup_name = request.form["workgroup"]
        workbook_name = request.form["workbook"]
    # validate user belongs to the workbook
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    abort_if_user_not_in_workbook(workgroup_name, workbook_name, workbook)


def edit_reaction(
    reaction: models.Reaction,
    request_source: str = "user",
    file_attachment: bool = False,
    file_object_for_deletion: models.ReactionDataFile = None,
):
    """
    Validates the active user has permission to edit the reaction. Protects against user edited HTML.
    Aborts process with a 401 error if validation is failed

    Args:
        reaction: the active reaction
        request_source: the origin of request. Either 'user' or 'server'
        file_attachment: Whether file attachments are being edited.
        file_object_for_deletion: If a file is being deleted, we check this file can be deleted
    """
    # validate user is in workbook
    workbook_persons = reaction.workbook.users
    workbook_users = [x.user for x in workbook_persons]
    if current_user not in workbook_users:
        abort(401)
    # validate the user is the creator
    if reaction.creator_person.user.email != current_user.email:
        abort(401)
    # validate the user is not deleting autogenerated files.
    if (
        file_attachment
        and request.method == "DELETE"
        and request_source == "user"
        and file_object_for_deletion.autogenerated is True
    ):
        abort(401)
    if file_attachment:
        return
    # validate the reaction is not locked, unless it is a file attachment being edited.
    if reaction.complete == "complete":
        abort(401)