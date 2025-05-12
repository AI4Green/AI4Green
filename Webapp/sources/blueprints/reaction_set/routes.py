from flask import render_template

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
