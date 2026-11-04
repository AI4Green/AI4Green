from flask import abort, jsonify
from flask_login import login_required
from sources.auxiliary import security_member_workgroup

from . import inventory_api_bp


@inventory_api_bp.route("/<workgroup_name>", methods=["GET"])
@login_required
def get_workgroup_inventory(workgroup_name):
    """
    Get the inventory for a workgroup: its storage locations and every compound
    it holds, each with its containers nested underneath.

    Currently returns mock data so the frontend can be built ahead of the real
    queries. The response shape follows the Location, WorkgroupCompound and
    InventoryItem models, so swapping in real data should not change the contract.
    """
    if not security_member_workgroup(workgroup_name):
        abort(404)

    # todo: replace with real queries once the inventory services exist
    return jsonify("{}")
