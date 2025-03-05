from flask import render_template, request
from flask_login import current_user, login_required

from . import retrosynthesis_bp


@retrosynthesis_bp.route("/retrosynthesis/", methods=["GET", "POST"])
@login_required
@retrosynthesis_bp.doc("sessionAuth")
def retrosynthesis():
    """
    Renders the retrosynthesis page

    Returns:
        flask.Response: renders the retrosynthesis template
    """
    current_user.retrosynthesis_smiles = request.form.get("smiles")
    return "", 204


@retrosynthesis_bp.route("/retrosynthesis_about", methods=["GET"])
@login_required
@retrosynthesis_bp.doc("sessionAuth")
def retrosynthesis_about():
    """
    Renders the retrosynthesis about page

    Returns:
        flask.Response: renders the retrosynthesis about template
    """
    return render_template("retrosynthesis_about.html")
