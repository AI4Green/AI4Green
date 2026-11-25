from apiflask import APIBlueprint

methods_api_bp = APIBlueprint("methods", __name__, url_prefix="/methods")

from . import routes, template_routes  # noqa: E402,F401
