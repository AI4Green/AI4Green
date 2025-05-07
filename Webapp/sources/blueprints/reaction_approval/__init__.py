from apiflask import APIBlueprint

reaction_approval_bp = APIBlueprint("reaction_approval", __name__)

from . import routes
