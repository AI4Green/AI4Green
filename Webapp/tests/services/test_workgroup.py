import pytest_mock
from flask import Flask
from flask.testing import FlaskClient
from sources import models
from sources.extensions import db
from tests.utils import login


def test_change_workgroup_name_success(
    app: Flask, client: FlaskClient, mocker: pytest_mock.MockerFixture
):
    """
    Tests the successful change of a workgroup's name when the new name verification matches the user input.
    """
    login(client)
    original_name = "Test-Workgroup"
    new_name = "New-Workgroup"
    mock_pi_access = mocker.patch("sources.decorators.services.workgroup.get_user_type")
    mock_pi_access.return_value = "principal_investigator"
    response = client.post(f"/change_workgroup_name/{original_name}/{new_name}")
    with app.app_context():
        new_name_wg = (
            db.session.query(models.WorkGroup)
            .filter(models.WorkGroup.name == new_name)
            .first()
        )
    assert new_name_wg is not None, "check the new workgroup name is not a null value"
    assert response.status_code == 200
    assert response.json == {
        "message": "Workgroup name successfully changed",
        "new_name": new_name,
    }, "The success message should match"


def test_change_workgroup_name_failure(client: FlaskClient):
    """
    Tests the failure response when the new name verification does not match the user input.
    """
    login(client)
    original_name = "Test-Workgroup"
    new_name = "****-Workgroup"
    error_message = "Name verification failed"
    response = client.post(f"/change_workgroup_name/{original_name}/{new_name}")
    db.session.commit.assert_not_called()
    assert response.status_code == 400
    assert response.json == {"error": error_message}, "The error message should match"
