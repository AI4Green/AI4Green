from flask import jsonify
from flask_login import current_user
from sources import db, models

from . import coshh_api_bp


@coshh_api_bp.route("/project_types", methods=["GET"])
def get_all_templates():
    return jsonify(
        [
            {"id": "uu", "name": "ggg", "stage": "Draft"},
            {"id": "tt", "name": "ddd", "stage": "Draft"},
        ]
    )
    # return (
    #     db.session.query(models.User)
    #     .filter(func.lower(models.User.email) == user_email.lower())
    #     .first()
    # )
