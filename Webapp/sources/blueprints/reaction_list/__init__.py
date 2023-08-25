from flask import Blueprint

reaction_list_bp = Blueprint("reaction_list", __name__)

from . import routes
