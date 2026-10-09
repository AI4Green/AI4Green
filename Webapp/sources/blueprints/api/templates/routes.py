from datetime import datetime
from tempfile import template

import pytz
import services.templates
from flask import jsonify, request
from flask_login import current_user
from sources import db, models

from . import templates_api_bp


@templates_api_bp.route("/", methods=["GET"])
def get_templates():
    template_type = request.args.get("template_type", None)
    print(template_type)
    template_list = []

    if template_type == "COSHH":
        template_list = services.templates.list_coshh()

    return jsonify(template_list)


@templates_api_bp.route("/<int:template_id>", methods=["GET"])
def get_template(template_id):
    template_object = models.Template.query.get(template_id)
    return jsonify(template_object.to_dict())


@templates_api_bp.route("/", methods=["POST"])
def create_new_template():
    data = request.get_json()
    source_id = data.get("source_id", None)
    name = data.get("name", None)
    description = data.get("description", None)
    template_type = data.get("templateType", None)

    new_template = {}

    # todo: handle errors if name or desc are missing
    # todo: include institution id per user

    # if no source id, create a blank template with default values
    if not source_id:
        if template_type == "COSHH":
            new_template = services.templates.create_new_coshh_template(
                name, description
            )

    return jsonify(new_template), 200


@templates_api_bp.route("/<int:template_id>/sections", methods=["GET"])
def get_template_sections(template_id):
    query = models.Template.query.get(template_id)
    return jsonify([x.to_dict() for x in query.sections])
