from flask_login import current_user

from . import users_api_bp


@users_api_bp.route("/me", methods=["GET", "POST"])
def me():
    # sparse route atm for passing login, needs to update permissions based on roles
    print(current_user.fullname)
    return {
        "fullName": current_user.fullname,
        "email": current_user.email,
        "permissions": [
            "CreateProjectTypes",
            "EditProjectTypes",
            "DeleteProjectTypes",
            "ViewProjectTypes",
        ],
    }
