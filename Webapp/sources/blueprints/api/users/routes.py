from flask_login import current_user
from sources.auxiliary import get_workgroups

from . import users_api_bp


@users_api_bp.route("/me", methods=["GET", "POST"])
def me():
    # sparse route atm for passing login, needs to update permissions based on roles
    return {
        "fullName": current_user.fullname,
        "userName": current_user.username,
        "email": current_user.email,
        "permissions": [
            "CreateProjectTypes",
            "EditProjectTypes",
            "DeleteProjectTypes",
            "ViewProjectTypes",
        ],
        "workgroups": get_workgroups(),
    }
