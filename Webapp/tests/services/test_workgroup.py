import pytest_mock
from flask import Flask
from flask.testing import FlaskClient
from sources import services
from sources.extensions import db
from tests.utils import add_workgroup, login


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
        new_name_wg = services.workgroup.from_name(new_name)
        assert new_name_wg is not None, "Check the workgroup is renamed successfully"
        assert response.status_code == 200
        assert response.json == {
            "message": "Workgroup name successfully changed",
            "new_name": new_name,
        }, "The success message should match"
        new_name_wg.name = original_name
        db.session.commit()


def test_change_workgroup_name_failure_invalid_symbols(app: Flask, client: FlaskClient):
    """
    Tests the failure response when the new name verification does not match the user input
    and contains invalid symbols.
    """
    login(client)
    original_name = "Test-Workgroup"
    new_name = "****-Workgroup"
    error_message = "Error: Name contains invalid symbols."
    response = client.post(f"/change_workgroup_name/{original_name}/{new_name}")
    with app.app_context():
        original_wg_name = services.workgroup.from_name(original_name)
    assert (
        original_wg_name is not None
    ), "Check the name change is unsuccessful due to an invalid symbols"
    assert response.status_code == 400
    assert response.json == {"error": error_message}, "The error message should match"


def test_change_workgroup_name_failure_duplicate(app: Flask, client: FlaskClient):
    """
    Tests the failure response when the new name verification does not match the user input
    and is a duplicate name.
    """
    login(client)
    with app.app_context():
        duplicate_workgroup = add_workgroup("duplicate")
        original_name = "Test-Workgroup"
        new_name = "duplicate"
        error_message = "Workgroup name already exists"
        response = client.post(f"/change_workgroup_name/{original_name}/{new_name}")
        original_wg_name = services.workgroup.from_name(original_name)
    assert (
        original_wg_name is not None
    ), "Check the new workgroup name change failed due to a duplicate existing"
    assert response.status_code == 400
    assert response.json == {"error": error_message}, "The error message should match"

    with app.app_context():
        db.session.delete(duplicate_workgroup)
        db.session.commit()
