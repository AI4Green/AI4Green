from flask import jsonify
from flask_login import login_required
from sources import services

from . import data_access_changes_bp


@data_access_changes_bp.route("/data_access_changes", methods=["GET"])
@login_required
def get_data_access_changes():
    person_id = services.person.from_current_user_email()
    changes = services.data_access_changes.get_changes_from_person(person_id)

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
