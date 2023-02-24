from flask import Blueprint

solvents_bp = Blueprint('solvents', __name__)

from sources.solvents import routes