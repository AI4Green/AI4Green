from apiflask import APIBlueprint

notifications_bp = APIBlueprint("notifications", __name__)

from . import routes
