from flask import abort
from flask_login import current_user
from sources import models


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
