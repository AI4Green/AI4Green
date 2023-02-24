from flask import Blueprint

novel_compound_bp = Blueprint('novel_compound', __name__)

from sources.novel_compound import routes
