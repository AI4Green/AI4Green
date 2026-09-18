from apiflask import APIBlueprint

sections_api_bp = APIBlueprint("sections", __name__, url_prefix="/sections")

from . import routes
