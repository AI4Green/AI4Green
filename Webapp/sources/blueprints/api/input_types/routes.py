from flask import jsonify
from sources import db, models

from . import input_types_api_bp


@input_types_api_bp.route("/", methods=["GET"])
def get_input_types():
    query = models.InputType.query.all()
    return jsonify([x.to_dict() for x in query])
