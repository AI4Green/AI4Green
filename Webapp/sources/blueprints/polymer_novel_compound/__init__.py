from flask import Blueprint

polymer_novel_compound_bp = Blueprint("polymer_novel_compound", __name__)

from . import routes
