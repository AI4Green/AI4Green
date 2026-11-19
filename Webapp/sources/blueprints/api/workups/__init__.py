from apiflask import APIBlueprint

workups_api_bp = APIBlueprint("workups", __name__, url_prefix="/workups")

from . import routes  # noqa: E402,F401
