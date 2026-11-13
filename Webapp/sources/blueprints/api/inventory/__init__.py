from apiflask import APIBlueprint

inventory_api_bp = APIBlueprint("inventory", __name__, url_prefix="/inventory")

from . import routes
