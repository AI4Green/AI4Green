from flask_login import current_user

from . import users_api_bp


@users_api_bp.route("/me", methods=["GET", "POST"])
def me():
    print(current_user.email)
    return {
        "message": "Connection Successful!",
        "user_status": "Authenticated",
        # "workgroup": workgroup
    }
