from apiflask import APIBlueprint

reset_password_bp = APIBlueprint("reset_password", __name__)

from . import routes
