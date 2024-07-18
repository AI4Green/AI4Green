from functools import wraps
from flask import redirect, url_for, flash, request
from flask_login import current_user, login_required
from sources import services
from sources.auxiliary import get_workgroups, get_workbooks


def _get_from_request(input_value, search_str):
    if search_str in request.args:
        return request.args.get(search_str)
    elif search_str in request.form.keys():
        return request.form.get(search_str)
    else:
        return input_value


def _is_demo():
    if request.form.get("demo") == "demo" or request.form.get("tutorial") == "tutorial":
        return True
    if request.args.get("demo") is None and request.args.get("tutorial") is None:
        return False
    return True


def principal_investigator_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        workgroup = _get_from_request(kwargs["workgroup"], "workgroup")
        if services.workgroup.get_user_type(workgroup, current_user) != "principal_investigator":
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        return f(*args, **kwargs)
    return decorated_function


def principal_investigator_or_senior_researcher_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        workgroup = _get_from_request(kwargs["workgroup"], "workgroup")
        if services.workgroup.get_user_type(workgroup, current_user) not in ["principal_investigator", "senior_researcher"]:
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        return f(*args, **kwargs)
    return decorated_function


def workgroup_member_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        workgroup = _get_from_request(kwargs["workgroup"], "workgroup")
        if workgroup not in get_workgroups():
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        return f(*args, **kwargs)
    return decorated_function


def workbook_member_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if not _is_demo():
            workgroup = _get_from_request(kwargs.get("workgroup"), "workgroup")
            workbook = _get_from_request(kwargs.get("workbook"), "workbook")
            if workbook not in get_workbooks(workgroup):
                flash("You do not have permission to view this page")
                return redirect((url_for("main.index")))
        return f(*args, **kwargs)
    return decorated_function

