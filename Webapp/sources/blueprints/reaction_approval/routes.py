from flask import (
    Response,
    current_app,
    flash,
    jsonify,
    redirect,
    render_template,
    request,
    url_for,
)
from flask_login import (  # protects a view function against anonymous users
    current_user,
    login_required,
)
from sources import models, services
from sources.auxiliary import (
    get_notification_number,
    get_workgroups,
    security_pi_workgroup,
)
from sources.decorators import principal_investigator_required, workbook_member_required

from . import reaction_approval_bp


@reaction_approval_bp.route("/new_request", methods=["POST"])
@login_required
def new_reaction_approval_request() -> Response:
    """
    Initiates a reaction approval request from reaction constructor

    Returns:
        flask.Response: success message
    """
    # CHECK FOR DUPLICATES AND PREVENT FOR MALICIOUS FRONT END CONTENT
    new_request = services.reaction.ReactionApprovalRequestSubmission()
    new_request.submit_request()

    return jsonify({"response": "success"})


@reaction_approval_bp.route("/resubmit_request", methods=["PATCH"])
@login_required
def resubmit_approval_request() -> Response:
    """
    Resubmits a reaction approval request after suggested changes are made

    Returns:
        flask.Response: success message
    """
    approval_request = models.ReactionApprovalRequest.query.get(
        request.json.get("approvalID")
    )
    current_request = services.reaction.ReactionApprovalRequestSubmission()
    current_request.resubmit_request(approval_request)

    return jsonify({"response": "success"})


@reaction_approval_bp.route("/request_response/<token>", methods=["GET", "POST"])
@login_required
def request_response(token: str) -> Response:
    """
    Render the reaction review interface if the user is authorized and the request is valid.

    If the approval request is still pending, the reaction sketcher interface is rendered for review.
    If the request has already been reviewed, the user is notified accordingly.

    Args:
        token (str): A secure token representing the reaction approval request.

    Returns:
        flask.Response: A rendered template for the reaction review, or a redirect with a flash message
        depending on authorization and approval request status.
    """
    (
        workgroup,
        workbook,
        reaction,
        reaction_approval_request,
    ) = services.email.verify_reaction_approval_request_token(token)
    addenda = services.reaction.get_addenda(reaction)

    # check user is a required approver AND a PI in the workgroup
    if (
        current_user.Person in reaction_approval_request.required_approvers
        and security_pi_workgroup(workgroup.name)
    ):
        if reaction_approval_request.status.value == "PENDING":
            return render_template(
                "sketcher_reload.html",
                reaction=reaction,
                load_status="loading",
                demo="not demo",
                active_workgroup=workgroup.name,
                active_workbook=workbook.name,
                tutorial="no",
                addenda=addenda,
                review=True,
            )
        else:
            flash("A review has already been submitted for this reaction!")
            return redirect(url_for("main.index"))

    else:
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))


@reaction_approval_bp.route("/approve", methods=["PATCH"])
@login_required
def review_reaction_approve() -> Response:
    """
    Handle a PATCH request to approve a reaction approval request and redirect the user to the home page

    Returns:
        flask.Response: A Flask JSON response containing a success message and a redirect URL.
    """
    approval_request = models.ReactionApprovalRequest.query.get(
        request.json.get("approvalID")
    )
    if approval_request is None:
        return jsonify(
            {
                "message": "This review request could not be found! Please contact the reaction creator.",
                "redirect_url": url_for("main.index"),
            }
        )

    request_status = services.reaction.ReactionApprovalRequestResponse(approval_request)
    request_status.approve()

    return jsonify(
        {"message": "Reaction approved!", "redirect_url": url_for("main.index")}
    )


@reaction_approval_bp.route("/reject", methods=["PATCH"])
@login_required
def review_reaction_reject() -> Response:
    """
    Handle a PATCH request to reject a reaction approval request and redirect the user to the home page

    Returns:
        flask.Response: A Flask JSON response containing a success message and a redirect URL.
    """
    approval_request = models.ReactionApprovalRequest.query.get(
        request.json.get("requestID")
    )
    if approval_request is None:
        return jsonify(
            {
                "message": "This review request could not be found! Please contact the reaction creator.",
                "redirect_url": url_for("main.index"),
            }
        )

    request_status = services.reaction.ReactionApprovalRequestResponse(approval_request)
    request_status.reject(request.json.get("comments"))

    return jsonify(
        {
            "message": "This reaction has been rejected. The creator will be notified!",
            "redirect_url": url_for("main.index"),
        }
    )


@reaction_approval_bp.route("/suggest_changes", methods=["PATCH"])
@login_required
def review_reaction_suggest_changes() -> Response:
    """
    Handle a PATCH request to suggest changes to reaction approval request and redirect the user to the home page

    Returns:
        flask.Response: A Flask JSON response containing a success message and a redirect URL.
    """
    approval_request = models.ReactionApprovalRequest.query.get(
        request.json.get("requestID")
    )
    if approval_request is None:
        return jsonify(
            {
                "message": "This review request could not be found! Please contact the reaction creator.",
                "redirect_url": url_for("main.index"),
            }
        )

    request_status = services.reaction.ReactionApprovalRequestResponse(approval_request)
    request_status.suggest_changes(request.json.get("comments"))

    return jsonify(
        {
            "message": "Your comments have been submitted! The reaction creator will be notified.",
            "redirect_url": url_for("main.index"),
        }
    )
