import json
from unittest.mock import MagicMock
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


def test_download_sends_a_zip_file(client: FlaskClient, mocker: MockerFixture):
    # Patch the MinIO client in the service module
    mock_client = mocker.patch("sources.services.audit_logs.client")

    # Create a mock object returned by list_objects
    mock_obj = MagicMock()
    mock_obj.object_name = "example.json"
    mock_client.list_objects.return_value = [mock_obj]

    # Create a mock file-like object for get_object
    mock_file = MagicMock()
    mock_content = json.dumps(
        {
            "person": 3,
            "workgroup": 2,
            "workbook": 2,
            "reaction": 4,
            "field_name": "Edited Reaction",
            "change_details": '{"reactant_physical_forms_text": {"new_value": ["Unknown", "-select-"], "old_value": ["Non-volatile liquid", "-select-"]}}',
            "date": "2025-06-11",
        }
    ).encode()
    mock_file.read.return_value = mock_content
    mock_client.get_object.return_value = mock_file

    login(client)
    expected_code = 200
    expected_mime_type = "application/zip"
    response = client.get(
        "/audit_log/Test-Workgroup/download",
        query_string={"topic": "test_topic"},
    )

    assert response.status_code == expected_code
    assert response.mimetype == expected_mime_type
