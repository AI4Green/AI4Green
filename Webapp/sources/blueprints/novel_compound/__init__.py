from apiflask import APIBlueprint

novel_compound_bp = APIBlueprint("novel_compound", __name__)

from . import routes
