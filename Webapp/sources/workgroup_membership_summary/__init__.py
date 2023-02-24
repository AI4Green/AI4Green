from flask import Blueprint

workgroup_membership_summary_bp = Blueprint('workgroup_membership_summary', __name__)

from sources.workgroup_membership_summary import routes
