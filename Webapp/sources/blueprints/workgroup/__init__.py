from apiflask import APIBlueprint

workgroup_bp = APIBlueprint("workgroup", __name__)

from . import routes
