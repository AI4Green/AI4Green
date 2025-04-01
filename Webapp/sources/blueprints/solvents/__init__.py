from apiflask import APIBlueprint

solvents_bp = APIBlueprint("solvents", __name__)

from . import routes
