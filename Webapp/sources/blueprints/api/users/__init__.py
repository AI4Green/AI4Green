from apiflask import APIBlueprint

users_api_bp = APIBlueprint("users", __name__, url_prefix="/users")

from . import routes
