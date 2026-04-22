from flask import jsonify
from sources import db, models

from . import sections_api_bp


@sections_api_bp.route("/<int:section_id>", methods=["GET"])
def get_template_sections(section_id):
    section = models.Section.query.get(section_id)
    return jsonify(section.to_dict)
