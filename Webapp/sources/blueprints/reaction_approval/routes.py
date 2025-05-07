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
    new_request = services.reaction.NewReactionApprovalRequest()
    new_request.submit_request()

    return jsonify({"response": "success"})


@reaction_approval_bp.route("/request_response/<token>", methods=["GET", "POST"])
@login_required
def reaction_approval_request_response(token: str) -> Response:
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
        flash("You do not have permission to view this page")
        return redirect(url_for("main.index"))


@reaction_approval_bp.route("/approve", methods=["PATCH"])
@login_required
@principal_investigator_required
def review_reaction_approve() -> Response:
    approval_request = models.ReactionApprovalRequest.query.get(
        request.json.get("approvalID")
    )
    request_status = services.reaction.ApprovalRequestStatus(approval_request)
    request_status.approve()

    flash("Reaction approved!")
    return jsonify(
        {"message": "Reaction approved!", "redirect_url": url_for("main.index")}
    )
