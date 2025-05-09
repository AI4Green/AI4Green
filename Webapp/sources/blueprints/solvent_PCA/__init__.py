from apiflask import APIBlueprint

solvent_PCA_bp = APIBlueprint("solvent_surfer", __name__)

from . import routes
