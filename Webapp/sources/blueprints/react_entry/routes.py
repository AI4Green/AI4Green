from flask import current_app, redirect
from flask_login import login_required

from . import react_entry_bp


@react_entry_bp.route("/coshh")
@login_required
def spa_coshh():
    react_url = current_app.config["REACT_APP_URL"]
    return redirect(react_url + "project-types/")
