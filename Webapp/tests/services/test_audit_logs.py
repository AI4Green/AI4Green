from flask import Flask
from flask.testing import FlaskClient

from tests.utils import login


def test_download_fails_without_topic(client: FlaskClient):
    login(client)
    expected_code = 400
    expected_msg = {"error": "topic is a required field"}
    response = client.get("/audit_log/download")

    assert response.status_code == expected_code
    assert isinstance(response.json, dict)
    assert "error" in response.json
    assert response.json["error"] == expected_msg["error"]


def test_download_fails_with_workgroup_as_non_integer(client: FlaskClient):
    login(client)
    expected_code = 400
    expected_msg = {"error": "workgroup must an integer"}
    response = client.get(
        "/audit_log/download",
        query_string={"topic": "test_topic", "workgroup": "hello"},
    )

    assert response.status_code == expected_code
    assert isinstance(response.json, dict)
    assert "error" in response.json
    assert response.json["error"] == expected_msg["error"]
