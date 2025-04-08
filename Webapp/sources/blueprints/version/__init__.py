from apiflask import APIBlueprint

version_bp = APIBlueprint("version", __name__)

from . import routes
