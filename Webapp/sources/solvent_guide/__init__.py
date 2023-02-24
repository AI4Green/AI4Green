from flask import Blueprint

solvent_guide_bp = Blueprint('solvent_guide', __name__)

from sources.solvent_guide import routes