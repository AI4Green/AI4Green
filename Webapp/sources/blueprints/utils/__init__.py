from flask import Blueprint

utils_bp = Blueprint("utils", __name__)

from . import routes
