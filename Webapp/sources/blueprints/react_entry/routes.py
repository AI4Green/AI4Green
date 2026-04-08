from flask import current_app, redirect
from flask_login import login_required

from . import react_entry_bp


@login_required
@react_entry_bp.route("/coshh", defaults={"path": ""})
def serve_react(path):
    react_url = current_app.config["REACT_APP_URL"]
    return redirect(react_url + "project-types/")
