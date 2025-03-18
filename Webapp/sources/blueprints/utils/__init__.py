from apiflask import APIBlueprint

utils_bp = APIBlueprint("utils", __name__)

from . import routes
