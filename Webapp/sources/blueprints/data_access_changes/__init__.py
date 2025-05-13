from apiflask import APIBlueprint

data_access_changes_bp = APIBlueprint("data_access_changes", __name__)

from . import routes
