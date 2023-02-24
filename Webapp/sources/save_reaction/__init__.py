from flask import Blueprint

save_reaction_bp = Blueprint('save_reaction', __name__)

from sources.save_reaction import routes
