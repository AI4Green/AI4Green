from apiflask import APIBlueprint

solvent_guide_bp = APIBlueprint("solvent_guide", __name__)

from . import routes
