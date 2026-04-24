from flask import jsonify, request
from flask_login import current_user
from sources import db, models

from . import sections_api_bp


@sections_api_bp.route("/", methods=["GET"])
def get_sections():
    # get sections per user
    query = (
        db.session.query(models.Section)
        .join(models.Template)
        .filter(models.Template.creator_id == current_user.id)
    )

    return jsonify([x.to_dict() for x in query.all()])


@sections_api_bp.route("/", methods=["POST"])
def save_new_section():
    data = request.get_json()

    name = data.get("name", None)
    sort_order = data.get("sortOrder", None)
    template_id = data.get("projectTypeId", None)

    new_section = models.Section.create(
        name=name,
        sort_order=sort_order,
        template_id=template_id,
        section_type_id=1,  # default for now, change later (possibly with modal?)
    )

    db.session.add(new_section)
    db.session.commit()

    return jsonify(new_section.to_dict())


@sections_api_bp.route("/<int:section_id>", methods=["GET"])
def get_section_by_id(section_id):
    section = models.Section.query.get(section_id)
    return jsonify(section.to_dict)


@sections_api_bp.route("/<int:section_id>/fields", methods=["GET"])
def get_section_fields(section_id):
    section = models.Section.query.get(section_id)
    return [x.to_dict() for x in section.fields]


@sections_api_bp.route("/<int:section_id>/fields", methods=["POST"])
def save_new_field(section_id):
    data = request.get_json()

    print(
        data
    )  # todo: if there is a field id, we need to update the field instead of creating a new one

    for f in data:
        print("F", f)
        # if an id exists, try to find it
        fid = f.get("id", None)
        if fid:
            query = models.Field.query.filter(models.Field.id == fid)
            query.update(
                {
                    "section_id": section_id,
                    "name": f.get("name"),
                    "sort_order": f.get("sortOrder"),
                    "input_type_id": f.get("inputType"),
                    "mandatory": f.get("mandatory"),
                }
            )

        else:
            new_field = models.Field.create(
                section_id=section_id,
                name=f.get("name"),
                sort_order=f.get("sortOrder"),
                input_type_id=f.get("inputType"),
                mandatory=f.get("mandatory"),
            )
            db.session.add(new_field)
    db.session.commit()

    return jsonify({"message": "success"}), 201
    # todo: handle multiple input fields, create new select_field_options?
