from flask import jsonify
from flask_login import login_required
from sources import db, models

from . import input_types_api_bp


@input_types_api_bp.route("/", methods=["GET"])
@login_required
def get_input_types():
    query = models.InputType.query.all()
    print([x.to_dict() for x in query])
    return jsonify([x.to_dict() for x in query])
