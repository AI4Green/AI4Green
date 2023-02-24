from flask import Blueprint

export_data_bp = Blueprint('export_data', __name__)

from sources.export_data import routes
