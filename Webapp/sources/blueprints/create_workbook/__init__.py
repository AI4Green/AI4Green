from apiflask import APIBlueprint

create_workbook_bp = APIBlueprint("create_workbook", __name__)

from . import routes
