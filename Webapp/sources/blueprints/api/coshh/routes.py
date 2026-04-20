from flask import jsonify
from flask_login import current_user
from sources import db, models

from . import coshh_api_bp


@coshh_api_bp.route("/project_types", methods=["GET"])
def get_all_templates():
    # templates = db.session.query(models.Template).filter()

    return jsonify(
        [
            {"id": "uu", "name": "ggg", "stage": "Draft", "permissions": "CanPublish"},
            {"id": "tt", "name": "ddd", "stage": "Draft", "permissions": "CanPublish"},
        ]
    )
    # return (
    #     db.session.query(models.User)
    #     .filter(func.lower(models.User.email) == user_email.lower())
    #     .first()
    # )
