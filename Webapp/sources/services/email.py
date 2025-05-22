from time import time
from typing import Dict, Optional, Tuple, Union

import jwt
from flask import current_app, render_template
from flask_login import current_user
from qrcode import QRCode
from sources import models, services
from sources.extensions import mail


def send_reaction_approval_request(
    user: models.User,
    reaction_approval_request: models.ReactionApprovalRequest,
    workgroup: models.WorkGroup,
    workbook: models.WorkBook,
    reaction: models.Reaction,
):
    """
    Send reaction approval request email to a given user. They have 70 days to respond
    Args:
        user (models.User): User to send to.
        reaction_approval_request (models.ReactionApprovalRequest): Reaction approval request.
        workgroup (models.WorkGroup): Workgroup the reaction belongs to
        workbook (models.WorkBook): Workbook the reaction belongs to
        reaction (models.Reaction): The reaction for approval
    """
    token = get_encoded_token(
        time_limit=60 * 60 * 24 * 70,  # 70 days
        arg_dict={
            "reaction_approval_request": reaction_approval_request.id,
            "workbook": workbook.id,
            "workgroup": workgroup.id,
            "reaction": reaction.id,
        },
    )
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Green Reaction Approval Request",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/reaction_approval_request.txt",
            reaction_approval_request=reaction_approval_request,
            user=user,
            workbook=workbook,
            workgroup=workgroup,
            token=token,
            protocol=protocol,
        ),
        html_body=render_template(
            "email/reaction_approval_request.html",
            reaction_approval_request=reaction_approval_request,
            user=user,
            workbook=workbook,
            workgroup=workgroup,
            token=token,
            protocol=protocol,
        ),
    )


def send_reaction_approval_response(
    user: models.User,
    reaction_approval_request: models.ReactionApprovalRequest,
):
    """
    Send reaction approval request response email to reaction creator. Updates the creator of a reaction on the verdict of the approval request
    Args:
        user (models.User): User to send to.
        reaction_approval_request (models.ReactionApprovalRequest): Reaction approval request.
    """

    mail.send_email(
        "Your Reaction Approval Request",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/reaction_approval_response.txt",
            reaction_approval_request=reaction_approval_request,
            user=user,
        ),
        html_body=render_template(
            "email/reaction_approval_response.html",
            reaction_approval_request=reaction_approval_request,
            user=user,
        ),
    )


def verify_reaction_approval_request_token(
    token: str,
) -> Optional[
    Dict[
        str,
        Union[
            models.WorkGroup,
            models.WorkBook,
            models.Reaction,
            models.ReactionApprovalRequest,
        ],
    ]
]:
    """
    Verify token link is valid and return workgroup, workbook and reaction ids for reaction to review

    Args:
        token: Token to verify.

    Returns:
        Tuple containing the workgroup, workbook, reaction and reaction_approval_request associated with the token.
        If token is invalid, returns None
    """
    decoded_token = decode_token(token)

    if not decoded_token:
        return None

    workgroup_id = decoded_token.get("workgroup")
    workbook_id = decoded_token.get("workbook")
    reaction_id = decoded_token.get("reaction")
    reaction_approval_request_id = decoded_token.get("reaction_approval_request")

    result = {
        "workgroup": models.WorkGroup.query.get(decoded_token.get(workgroup_id)),
        "workbook": models.WorkBook.query.get(decoded_token.get(workbook_id)),
        "reaction": models.Reaction.query.get(decoded_token.get(reaction_id)),
        "reaction_approval_request": models.ReactionApprovalRequest.query.get(
            decoded_token.get(reaction_approval_request_id)
        ),
    }

    # if any result fails return none
    if not all(result.values()):
        return None

    return result


def send_data_export_approval_request(
    user: models.User, data_export_request: models.DataExportRequest
):
    """
    Send data export request email to a given user. They have a week(or 604800 seconds) to respond
    Args:
        user: User to send to.
        data_export_request: request we are asking PI user to accept or deny
    """
    token = get_encoded_token(
        time_limit=604800,
        arg_dict={"data_export_request": data_export_request.id, "user": user.id},
    )
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Green Data Export Request",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/data_export_request.txt",
            data_export_request=data_export_request,
            user=user,
            token=token,
            protocol=protocol,
        ),
        html_body=render_template(
            "email/data_export_request.html",
            data_export_request=data_export_request,
            user=user,
            token=token,
            protocol=protocol,
        ),
    )


def send_data_export_ready_message(
    user: models.User, data_export_request: models.DataExportRequest
):
    """
    Send data export ready email message to a user. They have a week(or 604800 seconds) to respond

    Args:
        user: User to send to.
        data_export_request: The request which has been complete
    """
    token = get_encoded_token(
        time_limit=604800,
        arg_dict={"data_export_request": data_export_request.id, "user": user.id},
    )
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Green Data Export Ready",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/data_export_ready.txt",
            data_export_request=data_export_request,
            user=user,
            token=token,
            protocol=protocol,
        ),
        html_body=render_template(
            "email/data_export_ready.html",
            data_export_request=data_export_request,
            user=user,
            token=token,
            protocol=protocol,
        ),
    )


def get_encoded_token(time_limit: int, arg_dict: Dict) -> str:
    """
    Get token with expiry time.

    Args:
        time_limit: number of seconds before token expires. Use 600 for password resets and 7200 for email verification
        arg_dict: arguments to encode. Should contain either "password_reset" or "email_verification" key with user.id value

    Returns:
        A token string.
    """
    arg_dict["exp"] = time() + time_limit
    return jwt.encode(
        arg_dict,
        current_app.config["SECRET_KEY"],
        algorithm="HS256",
    )


def send_email_verification(user: models.User) -> None:
    """
    Send email verification to defined user.

    Args:
        user: User to send to.
    """
    token = get_encoded_token(time_limit=7200, arg_dict={"verify_email": user.id})
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Green Email Verification",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/email_verification.txt", user=user, token=token, protocol=protocol
        ),
        html_body=render_template(
            "email/email_verification.html", user=user, token=token, protocol=protocol
        ),
    )


def send_password_reset(user: models.User) -> None:
    """
    Send password reset email to a given user.

    Args:
        user: User to send to.
    """
    token = get_encoded_token(time_limit=600, arg_dict={"reset_password": user.id})
    protocol = get_protocol_type()
    mail.send_email(
        "AI4Green Reset Your Password",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[user.email],
        text_body=render_template(
            "email/reset_password_text.txt", user=user, token=token, protocol=protocol
        ),
        html_body=render_template(
            "email/reset_password_text.html", user=user, token=token, protocol=protocol
        ),
    )


def send_notification(person: models.Person):
    """
    Send notifications email to a user.

    Args:
        person: Person to send to.

    """
    protocol = get_protocol_type()
    mail.send_email(
        "You have a new AI4Green notification",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[person.user.email],
        text_body=render_template(
            "email/notification_text.txt", user=person.user, protocol=protocol
        ),
        html_body=render_template(
            "email/notification_text.html", user=person.user, protocol=protocol
        ),
    )


def verify_reset_password_token(token: str) -> Optional[models.User]:
    """
    Verify token link is valid and return user id.

    Args:
        token: Token to verify.

    Returns:
         User the token identifies.
    """
    decoded_token = decode_token(token)
    if decoded_token:
        user_id = decoded_token["reset_password"]
        return models.User.query.get(user_id)
    else:
        return None


def verify_data_export_token(
    token: str,
) -> Optional[Tuple[models.User, models.DataExportRequest]]:
    """
    Verify token link is valid and return user id.

    Args:
        token: Token to verify.

    Returns:
        Tuple containing User the token identifies and DataExportRequest associated with the token.
        If token is invalid, returns None for both User and DataExportRequest.
    """
    decoded_token = decode_token(token)

    if decoded_token:
        data_export_request_id = decoded_token.get("data_export_request")
        user_id = decoded_token.get("user")

        user = models.User.query.get(user_id)
        data_export_request = models.DataExportRequest.query.get(data_export_request_id)

        return user, data_export_request

    return None, None


def verify_qr_code_for_add_user_token(
    token: str,
) -> Optional[str]:
    """
    Verify token link is valid and return workgroup to add user to.

    Args:
        token: Token to verify.

    Returns:
        models.Workgroup, workgroup to add user to
    """
    decoded_token = decode_token(token)

    if decoded_token:
        workgroup_name = decoded_token.get("workgroup")
        return workgroup_name

    return None


def decode_token(token: str) -> Dict:
    try:
        decoded_token = jwt.decode(
            token, current_app.config["SECRET_KEY"], algorithms=["HS256"]
        )
        return decoded_token
    except jwt.exceptions.ExpiredSignatureError:
        return None


def send_password_reset_test(user: models.User) -> Tuple[str, str]:
    """
    Send password reset email in testing.

    Args:
        user: User to send to.

    Returns:
        A tuple of the rendered template, and the token.
    """
    token = get_encoded_token(time_limit=600, arg_dict={"reset_password": user.id})
    return (
        render_template("email/reset_password_text.html", user=user, token=token),
        token,
    )


def send_controlled_substance_alert(
    substance: str, location: Dict[str, str], reaction: models.Reaction
) -> str:
    """
    Send an email to AI4Green admin if use of a controlled substance has been detected in a country with
    a UK arms embargo.

    Args:
        substance: str, inchi of the substance that has been used
        location: Dict[str,str], the country where the substance has been used
        reaction: models.Reaction, the reaction that contains the controlled substance

    """
    compound = services.all_compounds.from_inchi(substance, reaction.workbook)
    mail.send_email(
        "Controlled Substance Alert",
        sender=current_app.config["MAIL_ADMIN_SENDER"],
        recipients=[
            current_app.config["EXPORT_CONTROL_EMAIL_ADDRESS"],
            current_app.config["MAIL_ADMIN_SENDER"],
        ],
        text_body=render_template(
            "email/controlled_substance_alert.txt",
            country=location["country"],
            substance_smiles=compound.smiles,
            email=current_user.email,
            city=location["city"],
            ip=location["IP_address"],
            substance_name=compound.name,
            substance_inchi=substance,
            reaction_id=reaction.reaction_id,
            reaction_name=reaction.name,
            reaction_smiles=reaction.reaction_smiles,
            workbook=reaction.workbook,
            date_created=reaction.time_of_creation,
            time_of_update=reaction.time_of_update,
        ),
        html_body=render_template(
            "email/controlled_substance_alert.html",
            country=location["country"],
            substance=substance,
            email=current_user.email,
            city=location["city"],
            ip=location["IP_address"],
            reaction_id=reaction.reaction_id,
            reaction_name=reaction.name,
            reaction_smiles=reaction.reaction_smiles,
            workbook=reaction.workbook,
            date_created=reaction.time_of_creation,
            time_of_update=reaction.time_of_update,
        ),
    )


# send reset email test
def send_notification_test(person: models.Person) -> str:
    """
    Send notification email in testing.

    Args:
        person: Person to send to.

    Returns:
         Rendered notification text template.
    """
    return render_template("email/notification_text.html", user=person.user)


def verify_encoded_token(token: str, identifier: str) -> models.User:
    """
    Verify token link is valid and return user id.

    Args:
        token: Token to verify.
        identifier: Dictionary key to query after token is decoded. Currently supports "reset_password" or "verify_email"

    Returns:
         User the token identifies.
    """
    try:
        user_id = jwt.decode(
            token, current_app.config["SECRET_KEY"], algorithms=["HS256"]
        )[identifier]
    except Exception:
        return
    return models.User.query.get(user_id)


def get_protocol_type() -> str:
    """
    The desired protocol type depends on the current deployment environment.
    http for local and https if the app is deployed on a remote server

    Returns:
          the active protocol type
    """
    return "http" if current_app.config["MAIL_USE_LOCAL"] == "local" else "https"
