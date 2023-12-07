from flask import Blueprint

solvent_PCA_bp = Blueprint("solvent_PCA", __name__)

from . import routes
