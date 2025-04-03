from apiflask import APIBlueprint

save_reaction_bp = APIBlueprint("save_reaction", __name__)

from . import routes
