from flask import Blueprint

workgroup_bp = Blueprint("workgroup", __name__)

from . import routes
