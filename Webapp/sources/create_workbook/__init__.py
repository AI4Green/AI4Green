from flask import Blueprint

create_workbook_bp = Blueprint('create_workbook', __name__)

from sources.create_workbook import routes
