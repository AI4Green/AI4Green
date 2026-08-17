from datetime import datetime

import pytz
from flask import jsonify, request
from flask_login import current_user
from sources import db, models

from . import coshh_api_bp


@coshh_api_bp.route("/templates/<int:template_id>/sections", methods=["GET"])
def get_template_sections(template_id):
    template = models.Template.query.get(template_id)
    return jsonify(template.sections)


# @coshh_api_bp.route("/templates/<int:template_id>/sections", methods=["GET"])
# def get_template_section_types(template_id):
#     print(template_id)
#
#     section_types = models.SectionType.query.all()
#     return jsonify([x.to_dict() for x in section_types])
