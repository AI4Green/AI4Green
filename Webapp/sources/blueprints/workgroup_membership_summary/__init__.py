from flask import Blueprint

workgroup_membership_summary_bp = Blueprint("workgroup_membership_summary", __name__)

from . import routes
