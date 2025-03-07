from apiflask import APIBlueprint

summary_bp = APIBlueprint("summary", __name__)

from . import routes
