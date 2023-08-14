from datetime import datetime

from flask import Response, request
from flask_login import login_required

from sources import models
from sources.extensions import db

from . import compound_data_error_report_bp


@compound_data_error_report_bp.route("/compound_data_error_report", methods=["POST"])
@login_required
def compound_data_error_report() -> Response:
    # must be logged in
    """Handles submitted compound data error reports by adding to the database"""
    compound_name = request.form["compoundName"]
    error_type = request.form["errorType"]
    additional_info = request.form["additionalInfo"]
    compound_id = request.form["compoundID"]
    compound = (
        db.session.query(models.Compound).filter(models.Compound == compound_id).first()
    )
    model = models.CompoundDataErrorReport(
        compound_name=compound_name,
        compound=compound,
        error_type=error_type,
        additional_info=additional_info,
        time=datetime.now(),
    )
    db.session.add(model)
    db.session.commit()
    return "report"
