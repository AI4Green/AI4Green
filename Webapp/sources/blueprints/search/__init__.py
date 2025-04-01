from apiflask import APIBlueprint

search_bp = APIBlueprint("search", __name__)

from . import routes
