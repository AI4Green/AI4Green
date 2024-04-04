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
