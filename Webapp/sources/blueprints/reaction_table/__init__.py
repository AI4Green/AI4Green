from apiflask import APIBlueprint

reaction_table_bp = APIBlueprint("reaction_table", __name__)

from . import routes
