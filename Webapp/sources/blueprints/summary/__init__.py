from flask import Blueprint

summary_bp = Blueprint("summary", __name__)

from . import routes
