import json
from unittest.mock import MagicMock
from flask import Flask
from flask.testing import FlaskClient
from pytest_mock import MockerFixture
from tests.utils import login


def test_download_fails_without_topic(client: FlaskClient):
    login(client)
    expected_code = 400
    expected_msg = {"error": "topic is a required field"}
    response = client.get("/audit_log/Test-Workgroup/download")

    assert response.status_code == expected_code
    assert isinstance(response.json, dict)
    assert "error" in response.json
    assert response.json["error"] == expected_msg["error"]


def test_download_fails_when_workgroup_non_existent(client: FlaskClient):
    login(client)
    expected_code = 302
    response = client.get(
        "/audit_log/Fake-Test-Workgroup/download",
        query_string={"topic": "test_topic"},
    )

    assert response.status_code == expected_code


def test_download_sends_a_zip_file(
    app: Flask, client: FlaskClient, mocker: MockerFixture
):
    with app.app_context():
        audit_message = {
            "full_name": "Gloria Testeban",
            "email": "test_user@test.com",
            "workgroup": 2,
            "workbook": 2,
            "reaction": 4,
            "field_name": "Edited Reaction",
            "change_details": '{"reactant_physical_forms_text": {"new_value": ["Unknown", "-select-"], "old_value": ["Non-volatile liquid", "-select-"]}}',
            "date": "2025-06-11",
        }
        # Mock the get_audit_logs service
        mock_service = mocker.patch("sources.services.audit_logs.get_audit_logs")
        mock_service.return_value = [audit_message]

        # Create a mock file-like object for get_object
        mock_file = MagicMock()
        mock_content = json.dumps(audit_message).encode()
        mock_file.read.return_value = mock_content

        login(client)
        expected_code = 200
        expected_mime_type = "application/zip"
        response = client.get(
            "/audit_log/Test-Workgroup/download",
            query_string={"topic": "reactionedithistory"},
        )

        assert response.status_code == expected_code
        assert response.mimetype == expected_mime_type


def test_get_human_readable_ids(app: Flask):
    logs = [
        {
            "person": 1,
            "workgroup": 1,
            "workbook": 1,
            "reaction": 4,
            "field_name": "Edited Reaction",
            "change_details": {
                "reactant_physical_forms_text": {
                    "new_value": ["Unknown", "-select-"],
                    "old_value": ["Non-volatile liquid", "-select-"],
                }
            },
            "date": "2025-06-11",
        }
    ]
    expected_output = [
        {
            "person": "Gloria Testeban",
            "workgroup": "Test-Workgroup",
            "workbook": "Test-Workbook",
            "reaction": 4,
            "field_name": "Edited Reaction",
            "change_details": {
                "reactant_physical_forms_text": {
                    "new_value": ["Unknown", "-select-"],
                    "old_value": ["Non-volatile liquid", "-select-"],
                }
            },
            "date": "2025-06-11",
        }
    ]

    with app.app_context():
        from sources.services.audit_logs import get_human_readable_ids

        get_human_readable_ids(logs=logs)
        # method alters logs in-place so it should now equal expected_logs
        assert logs == expected_output
