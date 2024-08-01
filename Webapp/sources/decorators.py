from typing import Callable, Any
from functools import wraps
from flask import redirect, url_for, flash, request, Response
from flask_login import current_user
from typing import Union
from sources import services
from sources.auxiliary import get_workgroups, get_workbooks


def _get_from_request(input_value: Union[str, None], search_str: str) -> Union[str, None]:
    """
    Used for decorators. Searches request.args or request.form for search_str. If that value exists in the request it is returned, else the
    input_value is returned
    Args:
        input_value: Value for variable from decorator function to return if request returns None
        search_str: Variable name to search in request (either workbook or workgroup)

    Returns:
        request.args.get(search_str) or request.form.get(search_str) if either is not None else input_value
    """
    if search_str in request.args:
        return request.args.get(search_str)
    elif search_str in request.form.keys():
        return request.form.get(search_str)
    else:
        return input_value


def _is_demo() -> bool:
    """
    checks whether demo or tutorial mode is active

    Returns:
        Bool: True if demo or tutorial is active else False

    """
    if request.form.get("demo") == "demo" or request.form.get("tutorial") == "tutorial":
        return True
    if request.args.get("demo") is None and request.args.get("tutorial") is None:
        return False
    return True


def principal_investigator_required(f: Callable[..., Response]) -> Callable[..., Response]:
    """
        decorates function f to check whether user is a Principal Investigator in the workgroup,
        and redirects to homepage if they are not.

        The workgroup to search is passed either through a URL variable (eg: url/<workgroup>/...) or extracted from a
        request using the _get_from_request function.

        NOTE: redirect may fail if f returns json to front end

        Returns:
            redirect to home if user is not PI or SR else f
        """
    @wraps(f)
    def decorated_function(*args: Any, **kwargs: Any):
        workgroup = _get_from_request(kwargs.get("workgroup"), "workgroup")
        if services.workgroup.get_user_type(workgroup, current_user) != "principal_investigator":
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        kwargs["workgroup"] = workgroup
        return f(*args, **kwargs)
    return decorated_function


def principal_investigator_or_senior_researcher_required(f: Callable[..., Response]) -> Callable[..., Response]:
    """
    Decorates function f to check whether user is a Principal Investigator or Senior Researcher in the workgroup,
    and redirects to homepage if they are not.

    The workgroup to search is passed either through a URL variable (eg: url/<workgroup>/...) or extracted from a
    request using the _get_from_request function.

    NOTE: redirect may fail if f returns json to front end

    Returns:
        redirect to home if user is not PI or SR else f
    """
    @wraps(f)
    def decorated_function(*args: Any, **kwargs: Any):
        workgroup = _get_from_request(kwargs.get("workgroup"), "workgroup")
        if services.workgroup.get_user_type(workgroup, current_user) not in ["principal_investigator", "senior_researcher"]:
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        kwargs["workgroup"] = workgroup
        return f(*args, **kwargs)
    return decorated_function


def workgroup_member_required(f: Callable[..., Response]) -> Callable[..., Response]:
    """
        Decorates function f to check whether user is member of the workgroup, and redirects to homepage if they are
        not.

        The workgroup to search is passed either through a URL variable (eg: url/<workgroup>/...) or extracted from a
        request using the _get_from_request function.

        NOTE: redirect may fail if f returns json to front end

        Returns:
            redirect to home if user is not PI or SR else f
        """
    @wraps(f)
    def decorated_function(*args: Any, **kwargs: Any):
        workgroup = _get_from_request(kwargs["workgroup"], "workgroup")
        if workgroup not in get_workgroups():
            flash("You do not have permission to view this page")
            return redirect(url_for("main.index"))
        kwargs["workgroup"] = workgroup
        return f(*args, **kwargs)
    return decorated_function


def workbook_member_required(f: Callable[..., Response]) -> Callable[..., Response]:
    """
        Decorates function f to check whether user is a member of the workbook,
        and redirects to homepage if they are not.

        The workgroup to search is passed either through a URL variable (eg: url/<workgroup>/...) or extracted from a
        request using the _get_from_request function.

        NOTE: redirect may fail if f returns json to front end

        Returns:
            redirect to home if user is not PI or SR else f
        """
    @wraps(f)
    def decorated_function(*args: Any, **kwargs: Any):
        if not _is_demo():
            workgroup = _get_from_request(kwargs.get("workgroup"), "workgroup")
            workbook = _get_from_request(kwargs.get("workbook"), "workbook")
            if workbook not in get_workbooks(workgroup):
                flash("You do not have permission to view this page")
                return redirect((url_for("main.index")))
            kwargs["workgroup"] = workgroup
            kwargs["workbook"] = workbook
        return f(*args, **kwargs)
    return decorated_function
