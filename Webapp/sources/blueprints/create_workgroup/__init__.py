from flask import Blueprint

create_workgroup_bp = Blueprint("create_workgroup", __name__)

from . import routes
