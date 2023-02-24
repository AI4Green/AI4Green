from flask import Blueprint

join_workgroup_bp = Blueprint('join_workgroup', __name__)

from sources.join_workgroup import routes
