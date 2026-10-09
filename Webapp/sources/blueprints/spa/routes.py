import os

from flask import current_app, render_template, send_from_directory
from sources import services

from . import spa_bp


def get_spa_dir():
    return os.path.join(
        current_app.root_path,
        "static",
        "spa",
    )


@spa_bp.route("/")
@spa_bp.route("/<path:path>")
def serve_spa_files(path=""):
    # for dev, use the vite dev server to serve react
    if current_app.config.get("VITE_DEV_SERVER"):
        vite_entry = services.spa.vite_entry
        return render_template("spa.html", vite_entry=vite_entry)

    # else use static build
    spa_dir = get_spa_dir()

    if path:
        full_path = os.path.join(spa_dir, path)

        if os.path.isfile(full_path):
            return send_from_directory(spa_dir, path)

    # React Router fallback
    return send_from_directory(spa_dir, "index.html")
