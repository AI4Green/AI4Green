from flask import Blueprint

reaction_table_bp = Blueprint("reaction_table", __name__)

from . import routes
