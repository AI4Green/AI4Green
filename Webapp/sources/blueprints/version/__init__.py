from flask import Blueprint

version_bp = Blueprint("version", __name__)

from . import routes
