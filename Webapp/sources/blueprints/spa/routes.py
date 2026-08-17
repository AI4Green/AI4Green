import os

from flask import current_app, send_from_directory

from . import spa_bp


def get_spa_dir():
    base_dir = os.path.abspath(os.path.join(current_app.root_path))
    return os.path.join(base_dir, "static", "spa")


@spa_bp.route("/")
def serve_spa_route():
    spa_dir = get_spa_dir()

    return send_from_directory(spa_dir, "index.html")


@spa_bp.route("/<path:path>")
def serve_spa_files(path):
    spa_dir = get_spa_dir()

    full_path = os.path.join(spa_dir, path)

    if os.path.exists(full_path):
        return send_from_directory(spa_dir, path)

    return send_from_directory(spa_dir, "index.html")
