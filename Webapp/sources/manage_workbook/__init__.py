from flask import Blueprint

manage_workbook_bp = Blueprint('manage_workbook', __name__)

from sources.manage_workbook import routes
