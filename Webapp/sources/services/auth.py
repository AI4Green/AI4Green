from flask import abort, request
from flask_login import current_user
from sources import models, services
from sources.auxiliary import abort_if_user_not_in_workbook


def reaction(permission_level: str, request_method: str, async_type: str = "AJAX"):
    """
    Authenticates user to either view or edit the reaction.
    permission_level: Takes value of 'edit' or 'view_only'
    request_method: Value of 'GET' changes behaviour, other strings all have same behaviour (e.g. POST, DELETE)
    """
    if permission_level == "view_only":
        view_files(request_method)
    if permission_level == "edit":
        reaction = services.reaction.get_current_from_request(async_type)
        edit_reaction(reaction, file_attachment=True)


def view_files(request_method):
    """Authenticates user as a workbook member or aborts. Gets the workgroup_name, workbook_name, and workbook."""
    if request_method == "GET":
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


def edit_reaction(reaction: models.Reaction, file_attachment=False):
    """
    In addition to frontend validation, backend validation protects against user edited HTML.
    Validates the active user is the creator and validates the reaction is not locked.
    Aborts process with a 400 error if validation is failed
    """
    # validate user is in workbook
    workbook_persons = reaction.workbook.users
    workbook_users = [x.user for x in workbook_persons]
    if current_user not in workbook_users:
        abort(401)
    # validate the user is the creator
    if reaction.creator_person.user.email != current_user.email:
        abort(401)
    if file_attachment:
        return
    # validate the reaction is not locked, unless it is a file attachment being edited.
    if reaction.complete == "complete":
        abort(401)
