from apiflask import APIBlueprint

search_api_bp = APIBlueprint("search", __name__, url_prefix="/search")

from . import routes
