from apiflask import APIBlueprint

audit_log_bp = APIBlueprint("audit_log", __name__)

from . import routes
