from functools import wraps
from flask import redirect, url_for, flash, request
from flask_login import current_user, login_required
from sources import services


def _get_workgroup():
    return request.args.get("workgroup")


def principal_investigator_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        workgroup = _get_workgroup()
        print(workgroup)
        if services.workgroup.get_user_type(workgroup) != "principal_investigator":
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        return f(*args, **kwargs)
    return decorated_function
