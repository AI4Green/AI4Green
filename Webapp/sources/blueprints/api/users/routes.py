from typing import Any, Dict

from flask_login import current_user, login_required
from sources import services
from sources.auxiliary import get_workgroups

from . import users_api_bp


@users_api_bp.route("/me", methods=["GET"])
@login_required
def me() -> Dict[str, Any]:
    """
    The current user's details and front-end permissions.

    COSHH templates are managed by H&S admins only. Workup templates are managed by
    workbook members, so they need no global permission here.
    """
    permissions = ["ViewCoshhTemplates"]
    if services.templates.is_hs_admin():
        permissions.extend(
            [
                "CreateCoshhTemplates",
                "EditCoshhTemplates",
                "DeleteCoshhTemplates",
            ]
        )
    # sparse route atm for passing login, needs to update permissions based on roles
    return {
        "fullName": current_user.fullname,
        "userName": current_user.username,
        "email": current_user.email,
        "permissions": permissions,
        "workgroups": get_workgroups(),
    }
