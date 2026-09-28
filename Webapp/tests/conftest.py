import json
import os
from unittest.mock import mock_open

import pytest
from sources import create_app, services
from sources.extensions import db
from tests import populate_database


@pytest.fixture(scope="session")
def app():
    """
    This fixture creates a temporary database and returns a Flask application
    with the database initialised.

    Returns:
        Flask application
    """
    app = create_app("test")

    with app.app_context():
        db.drop_all()
        db.create_all()
        # set up request context for location data to add users to db
        with app.test_request_context(
            json={
                "IP_address": "104.28.100.0",
                "country": "North Korea",
                "city": "Pyongyang",
            }
        ):
            populate_database.insert_test_data()

    yield app

    with app.app_context():
        db.drop_all()


@pytest.fixture
def client(app):
    """
    This fixture returns a Flask test client.

    Args:
        app: The Flask application
    """
    return app.test_client()


@pytest.fixture
def runner(app):
    """
    This fixture returns a Flask test CLI runner.

    Args:
        app: The Flask application
    """
    return app.test_cli_runner()


@pytest.fixture(autouse=True)
def mock_vite_manifest(monkeypatch):
    manifest = {
        "src/main.jsx": {
            "file": "assets/main-test.js",
            "css": ["assets/main-test.css"],
            "isEntry": True,
        },
        "src/entries/coshh-form.jsx": {
            "file": "assets/coshh-form-test.js",
            "css": ["assets/coshh-form-test.css"],
            "isEntry": True,
        },
    }

    real_open = open

    def mocked_open(file, *args, **kwargs):
        if (
            os.fspath(file)
            .replace("\\", "/")
            .endswith("/static/spa/.vite/manifest.json")
        ):
            return mock_open(read_data=json.dumps(manifest))()

        return real_open(file, *args, **kwargs)

    monkeypatch.setattr(
        "builtins.open",
        mocked_open,
    )
