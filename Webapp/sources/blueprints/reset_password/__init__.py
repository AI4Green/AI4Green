from flask import Blueprint

reset_password_bp = Blueprint("reset_password", __name__)

from . import routes
