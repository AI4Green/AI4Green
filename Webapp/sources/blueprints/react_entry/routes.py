from flask import current_app, redirect
from flask_login import current_user, login_required

from . import react_entry_bp


@react_entry_bp.route("/coshh")
@login_required
def spa_coshh():
    print(current_user.email)
    react_url = current_app.config["REACT_APP_URL"]
    return redirect(react_url + "project-types/")
