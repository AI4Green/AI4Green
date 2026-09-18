from apiflask import APIBlueprint

input_types_api_bp = APIBlueprint("input_types", __name__, url_prefix="/input-types")

from . import routes
