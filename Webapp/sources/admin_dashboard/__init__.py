from flask import Blueprint

admin_dashboard_bp = Blueprint('admin_dashboard', __name__)

from sources.admin_dashboard import routes