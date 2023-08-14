from flask import Blueprint

delete_profile_bp = Blueprint("delete_profile", __name__)

from . import routes
