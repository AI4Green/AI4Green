from apiflask import APIBlueprint

reaction_list_bp = APIBlueprint("reaction_list", __name__)

from . import routes
