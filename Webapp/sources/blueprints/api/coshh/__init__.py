from apiflask import APIBlueprint

coshh_api_bp = APIBlueprint("coshh", __name__, url_prefix="/coshh")

from . import routes
