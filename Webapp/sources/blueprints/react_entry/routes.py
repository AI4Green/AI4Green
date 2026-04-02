from flask import render_template

from . import spa_bp


@spa_bp.route("/", defaults={"path": ""})
@spa_bp.route("/<path:path>")
def serve_react(path):
    # this renders the same HTML any /spa/<> route
    # allows React Router to take over internal navigation.
    # i think
    return render_template("react_entry.html")
