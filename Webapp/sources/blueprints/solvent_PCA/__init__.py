from flask import Blueprint

solvent_PCA_bp = Blueprint("solvent_surfer", __name__)

from . import routes
