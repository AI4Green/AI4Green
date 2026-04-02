from flask import render_template

from . import react_entry_bp


@react_entry_bp.route("/", defaults={"path": ""})
@react_entry_bp.route("/<path:path>")
def serve_react(path):
    # this renders the same HTML any /spa/<> route
    # allows React Router to take over internal navigation.
    # i think
    return render_template("react_entry.html")
