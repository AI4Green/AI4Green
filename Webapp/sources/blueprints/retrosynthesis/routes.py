from flask import request
from flask_login import current_user, login_required

from . import retrosynthesis_bp


@retrosynthesis_bp.route("/retrosynthesis/", methods=["GET", "POST"])
@login_required
def retrosynthesis():
    current_user.retrosynthesis_smiles = request.form.get("smiles")
    return "", 204
