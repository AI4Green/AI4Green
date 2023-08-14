from flask import Blueprint

standard_landing_page_bp = Blueprint("standard_landing_page", __name__)

from . import routes
