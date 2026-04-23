from apiflask import APIBlueprint

fields_api_bp = APIBlueprint("fields", __name__, url_prefix="/fields")

from . import routes
