from typing import Dict

from flask.testing import FlaskClient
from sources import models
from sources.extensions import db


def login_response(client: FlaskClient, username: str, password: str):
    """Function to make login call to test client"""
    return client.post(
        "/home", data={"username": username, "password": password, "submit": True}
    )


def login(
    client: FlaskClient, username: str = "test_username", password: str = "test_pw"
):
    """Frequently used function to log a user in to gain access to test @login_required functions."""

    response = login_response(client, username, password)
    assert response.status_code == 302 and response.location == "/home", "login failed"

    # Ensure that the user is logged in by checking if the session contains user information
    with client.session_transaction() as session:
        assert "_user_id" in session, "user should be logged in"


def assert_expected_values(actual_data: Dict, expected_data: Dict):
    for key, expected_value in expected_data.items():
        assert (
            actual_data[key] == expected_value
        ), f"{expected_value} not found in data for key: {key}"


def add_workgroup(name: str) -> models.WorkGroup:
    """Adds a workgroup with default parameters to the database.
    Args:
        name: str, the name of the workgroup to be created.
    Returns:
        workgroup: WorkGroup object, the created workgroup."""

    p1 = db.session.query(models.Person).first()
    institution1 = db.session.query(models.Institution).first()
    workgroup = models.WorkGroup.create(
        name=name, institution=institution1.id, principal_investigator=[p1]
    )
    db.session.commit()
    return workgroup
