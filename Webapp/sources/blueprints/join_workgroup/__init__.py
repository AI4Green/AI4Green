from apiflask import APIBlueprint

join_workgroup_bp = APIBlueprint("join_workgroup", __name__)

from . import routes
