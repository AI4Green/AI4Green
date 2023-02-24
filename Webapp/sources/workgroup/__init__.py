from flask import Blueprint

workgroup_bp = Blueprint('workgroup', __name__)

from sources.workgroup import routes