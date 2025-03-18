from apiflask import APIBlueprint

admin_dashboard_bp = APIBlueprint("admin_dashboard", __name__)

from . import routes
