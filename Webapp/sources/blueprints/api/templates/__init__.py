from apiflask import APIBlueprint

templates_api_bp = APIBlueprint("templates", __name__, url_prefix="/templates")

from . import routes
