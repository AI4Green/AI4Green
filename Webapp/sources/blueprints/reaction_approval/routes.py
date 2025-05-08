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
def new_reaction_approval_request() -> Response:
    """
    Initiates a reaction approval request from reaction constructor
    """
    # CHECK FOR DUPLICATES AND PREVENT FOR MALICIOUS FRONT END CONTENT
    new_request = services.reaction.ReactionApprovalRequestSubmission()
    new_request.submit_request()

    return jsonify({"response": "success"})


@reaction_approval_bp.route("/resubmit_request", methods=["PATCH"])
def resubmit_approval_request() -> Response:
    """
    Resubmits a reaction approval request after suggested changes are made
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
# @principal_investigator_required
def review_reaction_approve() -> Response:
    approval_request = models.ReactionApprovalRequest.query.get(
        request.json.get("approvalID")
    )
    request_status = services.reaction.ReactionApprovalRequestResponse(approval_request)
    request_status.approve()

    return jsonify(
        {"message": "Reaction approved!", "redirect_url": url_for("main.index")}
    )


@reaction_approval_bp.route("/reject", methods=["PATCH"])
@login_required
# @principal_investigator_required
def review_reaction_reject() -> Response:
    approval_request = models.ReactionApprovalRequest.query.get(
        request.json.get("requestID")
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
    approval_request = models.ReactionApprovalRequest.query.get(
        request.json.get("requestID")
    )
    request_status = services.reaction.ReactionApprovalRequestResponse(approval_request)
    request_status.suggest_changes(request.json.get("comments"))

    return jsonify(
        {
            "message": "Your comments have been submitted! The reaction creator will be notified.",
            "redirect_url": url_for("main.index"),
        }
    )
