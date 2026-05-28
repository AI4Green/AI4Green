from apiflask import APIBlueprint

spa_bp = APIBlueprint("spa", __name__, url_prefix="/spa")

from . import routes
