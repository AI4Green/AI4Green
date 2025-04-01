from apiflask import APIBlueprint

delete_profile_bp = APIBlueprint("delete_profile", __name__)

from . import routes
