from apiflask import APIBlueprint

data_access_history_bp = APIBlueprint("data_access_history", __name__)
data_export_history_bp = APIBlueprint("data_export_history", __name__)

from . import routes
