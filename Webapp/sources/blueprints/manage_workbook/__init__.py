from apiflask import APIBlueprint

manage_workbook_bp = APIBlueprint("manage_workbook", __name__)

from . import routes
