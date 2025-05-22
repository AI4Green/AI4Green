from flask import jsonify, request
from flask_login import login_required
from sources import services
from sources.decorators import principal_investigator_required

from . import (
    data_access_history_bp,
    data_access_history_workgroup_bp,
    data_export_history_bp,
    data_export_history_workgroup_bp,
)


@data_access_history_bp.route("/data_access_history", methods=["GET"])
@login_required
def get_data_access_history():
    person = services.person.from_current_user_email()
    changes = services.data_access_history.get_history_from_person(person.id)

    result = [
        {
            "workgroup": services.workgroup.from_id(
                change.workgroup_id
            ).name,  # convert id to name
            "workbook": services.workbook.get(change.workbook_id).name
            if change.workbook_id
            else None,
            "old_role": change.old_role,
            "new_role": change.new_role,
            "time": change.time,
        }
        for change in changes
    ]

    return jsonify(result)


@data_access_history_workgroup_bp.route(
    "/data_access_history_workgroup", methods=["GET"]
)
@login_required
@principal_investigator_required
def get_data_access_history_workgroup():
    workgroup_name = request.args.get("workgroupName")
    workgroup_id = services.workgroup.from_name(workgroup_name).id
    changes = services.data_access_history.get_history_from_workgroup(workgroup_id)

    result = [
        {
            "person": change.person_id,
            "workbook": services.workbook.get(
                change.workbook_id
            ).name  # convert id to name
            if change.workbook_id
            else None,
            "old_role": change.old_role,
            "new_role": change.new_role,
            "time": change.time,
        }
        for change in changes
    ]

    return jsonify(result)


@data_export_history_bp.route("/data_export_history", methods=["GET"])
@login_required
def get_data_export_history():
    person = services.person.from_current_user_email()
    approved_exports = services.data_export.requests.get_approved_exports_from_person(
        person.id
    )

    result = [
        {
            "workgroup": services.workgroup.from_id(export.workgroup).name,
            "workbook": services.data_export.requests.get_workbook_name_from_id(
                export.id
            ),
            "time": export.time_of_request,
            "reactions": services.data_export.requests.get_reaction_id_list_from_id(
                export.id
            ),
        }
        for export in approved_exports
    ]

    return jsonify(result)


@data_export_history_workgroup_bp.route(
    "/data_export_history_workgroup", methods=["GET"]
)
@login_required
@principal_investigator_required
def get_data_export_history_workgroup():
    workgroup_name = request.args.get("workgroupName")
    workgroup_id = services.workgroup.from_name(workgroup_name).id
    approved_exports = (
        services.data_export.requests.get_approved_exports_from_workgroup(workgroup_id)
    )

    result = [
        {
            "person": export.requestor,
            "workbook": services.data_export.requests.get_workbook_name_from_id(
                export.id  # convert id to name
            ),
            "time": export.time_of_request,
            "reactions": services.data_export.requests.get_reaction_id_list_from_id(
                export.id
            ),
        }
        for export in approved_exports
    ]

    return jsonify(result)
