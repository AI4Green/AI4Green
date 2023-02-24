from flask import Blueprint

manage_workgroup_bp = Blueprint('manage_workgroup', __name__)

from sources.manage_workgroup import routes
