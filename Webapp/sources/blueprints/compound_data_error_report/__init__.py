from apiflask import APIBlueprint

compound_data_error_report_bp = APIBlueprint("compound_data_error_report", __name__)

from . import routes
