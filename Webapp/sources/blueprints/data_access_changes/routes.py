from flask import jsonify
from flask_login import login_required
from sources import services

from . import data_access_changes_bp, data_export_changes_bp


@data_access_changes_bp.route("/data_access_changes", methods=["GET"])
@login_required
def get_data_access_changes():
    person = services.person.from_current_user_email()
    changes = services.data_access_changes.get_changes_from_person(person.id)

    result = [
        {
            "workgroup": services.workgroup.from_id(
                change.workgroup_id
            ).name,  # convert id to name
            "workbook": services.workbook.get(change.workbook_id).name,
            "old_role": change.old_role,
            "new_role": change.new_role,
            "time": change.time,
        }
        for change in changes
    ]

    return jsonify(result)


@data_export_changes_bp.route("/data_export_changes", methods=["GET"])
@login_required
def get_data_export_changes():
    person = services.person.from_current_user_email()
    approved_exports = services.data_export.requests.get_approved_exports_from_person(
        person.id
    )

    result = [
        {
            "workgroup": services.data_export.requests.get_workgroup_name_from_id(
                export.data_export_request_id
            ),  # convert id to name
            "workbook": services.data_export.requests.get_workbook_name_from_id(
                export.data_export_request_id
            ),
            "time": services.data_export.requests.get_export_from_id(
                export.data_export_request_id
            ).time_of_request,
            "reactions": services.data_export.requests.get_reaction_id_list_from_id(
                export.data_export_request_id
            ),
        }
        for export in approved_exports
    ]

    return jsonify(result)
