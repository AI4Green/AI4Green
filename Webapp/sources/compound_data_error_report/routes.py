from flask import request
from flask_login import login_required
from sources.compound_data_error_report import compound_data_error_report_bp
from sources import db
from datetime import datetime
from pony.orm import select


@compound_data_error_report_bp.route('/compound_data_error_report', methods=['POST'])
@login_required
def compound_data_error_report():
	# must be logged in
	"""Handles submitted compound data error reports by adding to the database"""
	compound_name = request.form['compoundName']
	error_type = request.form['errorType']
	additional_info = request.form['additionalInfo']
	compound_id = request.form['compoundID']
	compound = select(x for x in db.Compound if x.id == compound_id).first()
	db.CompoundDataErrorReport(compound_name=compound_name, compound=compound, error_type=error_type, additional_info=additional_info, time=datetime.now())
	return 'report'

