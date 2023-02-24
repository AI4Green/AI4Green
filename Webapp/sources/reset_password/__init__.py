from flask import Blueprint

reset_password_bp = Blueprint('reset_password', __name__)

from sources.reset_password import routes