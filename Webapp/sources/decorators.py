from functools import wraps
from flask import redirect, url_for, flash, request
from flask_login import current_user, login_required
from sources import services
from sources.auxiliary import get_workgroups


def _get_workgroup(workgroup):
    if "workgroup" in request.args:
        return request.args.get("workgroup")
    else:
        return workgroup


def _get_workbook(workbook):
    if "workbook" in request.args:
        return request.args.get("workbook")
    else:
        return workbook


def principal_investigator_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        workgroup = _get_workgroup(kwargs["workgroup"])
        if services.workgroup.get_user_type(workgroup, current_user) != "principal_investigator":
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        return f(*args, **kwargs)
    return decorated_function


def workgroup_member_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        workgroup = _get_workgroup(kwargs["workgroup"])
        if workgroup not in get_workgroups():
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        return f(*args, **kwargs)
    return decorated_function


def workbook_member_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):

        workbook = _get_workbook(kwargs["workbook"])
        if workbook not in get_workgroups():
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        return f(*args, **kwargs)
    return decorated_function
