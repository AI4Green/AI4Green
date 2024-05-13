from flask import Blueprint

email_verification_bp = Blueprint("email_verification", __name__)

from . import routes
