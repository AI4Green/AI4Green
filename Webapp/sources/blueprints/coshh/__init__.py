from apiflask import APIBlueprint

coshh_bp = APIBlueprint("coshh", __name__)

from . import routes
