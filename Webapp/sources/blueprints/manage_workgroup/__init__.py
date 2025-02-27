from apiflask import APIBlueprint

manage_workgroup_bp = APIBlueprint("manage_workgroup", __name__)

from . import routes
