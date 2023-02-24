from flask import Blueprint

compound_data_error_report_bp = Blueprint('compound_data_error_report', __name__)

from sources.compound_data_error_report import routes