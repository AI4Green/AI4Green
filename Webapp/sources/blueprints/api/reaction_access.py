"""
Access checks shared by the API blueprints for the cards on the reaction page
(workups, methods).

Reactions are identified by the names the page uses: workgroup, workbook and
reaction code. Anyone who isn't a member of the workbook gets a 404, so the routes
don't reveal which workbooks or reactions exist. Only a reaction's creator can
change it, and only until it is locked (403 otherwise).
"""

from apiflask import abort
from sources import models, services


def member_workbook(workgroup_name: str, workbook_name: str) -> models.WorkBook:
    """
    Gets a workbook the current user is a member of, or stops with a 404.

    Args:
        workgroup_name: The workgroup name from the URL.
        workbook_name: The workbook name from the URL.

    Returns:
        The workbook.
    """
    workbook = services.workbook.get_workbook_from_group_book_name_combination(
        workgroup_name, workbook_name
    )
    if workbook is None or not services.workbook.is_member(workbook):
        abort(404)
    return workbook


def member_reaction(
    workgroup_name: str, workbook_name: str, reaction_id: str
) -> models.Reaction:
    """
    Gets a reaction in a workbook the current user is a member of, or stops with
    a 404.

    Args:
        workgroup_name: The workgroup name from the URL.
        workbook_name: The workbook name from the URL.
        reaction_id: The reaction code from the URL, e.g. "WB1-001".

    Returns:
        The reaction.
    """
    reaction = services.reaction.get_from_names(
        workgroup_name, workbook_name, reaction_id
    )
    if reaction is None:
        abort(404)
    return reaction


def can_edit(reaction: models.Reaction) -> bool:
    """
    Whether the current user may change the reaction: its creator, while the
    reaction isn't locked (the same rule as the reaction page's autosave).

    Args:
        reaction: The reaction.

    Returns:
        True if the user may change it.
    """
    return (
        services.reaction.check_current_user_is_creator(reaction)
        and reaction.complete != "complete"
    )


def editable_reaction(
    workgroup_name: str, workbook_name: str, reaction_id: str, locked_message: str
) -> models.Reaction:
    """
    Gets a reaction the current user may change, or stops with a 404 (not a
    member) or 403 (not the creator, or locked).

    Args:
        workgroup_name: The workgroup name from the URL.
        workbook_name: The workbook name from the URL.
        reaction_id: The reaction code from the URL.
        locked_message: The 403 message.

    Returns:
        The reaction.
    """
    reaction = member_reaction(workgroup_name, workbook_name, reaction_id)
    if not can_edit(reaction):
        abort(403, locked_message)
    return reaction
