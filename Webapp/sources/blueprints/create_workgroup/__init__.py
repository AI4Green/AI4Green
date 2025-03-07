from apiflask import APIBlueprint

create_workgroup_bp = APIBlueprint("create_workgroup", __name__)

from . import routes
