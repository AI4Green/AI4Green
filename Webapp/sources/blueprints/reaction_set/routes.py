import requests
from flask import current_app, render_template, request

from . import reaction_set_bp


@reaction_set_bp.route("/reaction_set")
def reaction_set():
    # to do
    # add workbook, workgroup to reaction set page
    # fix atuosave sketcher to only update reaction table in set mode (do we want to try autosaving?)
    # then fix apply to well and apply to all
    # colours for unsaved/edited wells
    # import from reactwise (try at least one small bit of data)
    #
    return render_template("reaction_set.html")


@reaction_set_bp.route("/click_and_drag")
def click_and_drag():
    return render_template("click-and-drag.html")


@reaction_set_bp.route("/import_from_reactwise", methods=["GET", "POST"])
def import_from_reactwise():
    step_id = request.json["step_id"]

    url = f"https://api.reactwise.com/views/step-experiments/{step_id}"

    payload = {}
    headers = {"Authorization": "Bearer " + current_app.config["REACTWISE_API_KEY"]}

    response = requests.request("GET", url, headers=headers, data=payload)

    print(response.text)
