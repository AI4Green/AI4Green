from flask import Blueprint

solvent_guide_bp = Blueprint("solvent_guide", __name__)

from . import routes
