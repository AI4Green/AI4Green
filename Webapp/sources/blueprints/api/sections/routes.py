from flask import jsonify, request
from flask_login import current_user
from sources import db, models

from . import sections_api_bp


@sections_api_bp.route("/", methods=["GET"])
def get_sections():
    template_id = request.args.get("template_id", None)

    # get sections per user
    query = db.session.query(models.Section).filter(
        models.Template.creator_id == current_user.id
    )

    if template_id:
        query.filter(models.Section.template_id == template_id)

    return jsonify([x.to_dict() for x in query.all()])


@sections_api_bp.route("/<int:section_id>", methods=["GET"])
def get_section_by_id(section_id):
    section = models.Section.query.get(section_id)
    return jsonify(section.to_dict)
