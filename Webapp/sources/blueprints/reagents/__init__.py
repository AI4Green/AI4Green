from apiflask import APIBlueprint

reagents_bp = APIBlueprint("reagents", __name__)

from . import routes
