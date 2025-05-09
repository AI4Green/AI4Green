from apiflask import APIBlueprint

retrosynthesis_bp = APIBlueprint("retrosynthesis", __name__)

from . import routes
