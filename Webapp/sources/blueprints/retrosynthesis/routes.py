from flask import render_template, request
from flask_login import current_user, login_required

from . import retrosynthesis_bp


@retrosynthesis_bp.route("/retrosynthesis/", methods=["GET", "POST"])
@login_required
def retrosynthesis():
    current_user.retrosynthesis_smiles = request.form.get("smiles")
    return "", 204


@retrosynthesis_bp.route("/retrosynthesis_about", methods=["GET"])
@login_required
def retrosynthesis_about():
    return render_template("retrosynthesis_about.html")
