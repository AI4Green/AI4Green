from flask import Blueprint

reagents_bp = Blueprint("reagents", __name__)

from . import routes
