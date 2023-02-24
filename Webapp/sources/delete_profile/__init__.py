from flask import Blueprint

delete_profile_bp = Blueprint('delete_profile', __name__)

from sources.delete_profile import routes
