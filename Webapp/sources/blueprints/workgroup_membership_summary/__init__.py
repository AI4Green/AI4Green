from apiflask import APIBlueprint

workgroup_membership_summary_bp = APIBlueprint("workgroup_membership_summary", __name__)

from . import routes
