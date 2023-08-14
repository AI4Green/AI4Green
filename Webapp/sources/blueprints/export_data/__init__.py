from flask import Blueprint

export_data_bp = Blueprint("export_data", __name__)

from . import routes
