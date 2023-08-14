from flask import Blueprint

manage_workbook_bp = Blueprint("manage_workbook", __name__)

from . import routes
