from flask import Response, jsonify
from flask_login import login_required
from sources import services

from . import utils_bp


@utils_bp.route("/updated_workgroup_dropdown", methods=["POST"])
@login_required
@utils_bp.doc(security="sessionAuth")
def updated_workgroup_dropdown() -> Response:
    """
    When the workgroup dropdown is updated, gets all workbooks which belong to the selected workgroup

    Returns:
        flask.Response: returns the workbooks belonging to the workgroup as a json object

    """
    # institution = request.form['institution']
    workbooks = services.interactive.update_workbook_dropdown()
    return jsonify({"workbooks": workbooks})
