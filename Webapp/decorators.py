from functools import wraps
from flask import redirect, url_for, flash
from sources import services

def workgroup_principal_investigator_required(workgroup):
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            if services.workgroup.get_user_type(workgroup.name) != "principal_investigator":
                flash("You do not have permission to view this page")
                return redirect(url_for("main.index"))
