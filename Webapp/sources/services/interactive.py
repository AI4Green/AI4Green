from flask import request
from sources.auxiliary import get_workbooks


def update_workbook_dropdown():
    """
    When workgroup changes, the workbook dropdown also should change. Some pages support all. Others don't.
    """
    workgroup = request.form["workgroup"]
    workbooks = get_workbooks(workgroup)
    if request.form.get("origin_page") != "export_data":
        workbooks.insert(0, "All")
    return workbooks