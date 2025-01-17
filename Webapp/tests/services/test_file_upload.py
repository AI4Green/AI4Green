import io
import os

import pytest
import pytest_mock
from flask.testing import FlaskClient
from tests.utils import login

IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.mark.skipif(IN_GITHUB_ACTIONS, reason="Azure blob storage is not setup in CI")
def test_reaction_file_attachment_upload_successful(
    client: FlaskClient, mocker: pytest_mock.MockerFixture
):
    """Tests the response when we try and upload a file attachment to a reaction"""
    login(client)
    # Define form data to send in the request
    form_data = {
        "workgroup": "Test-Workgroup",
        "workbook": "Test-Workbook",
        "reactionID": "TW1-001",
        "file": (io.BytesIO(b"abcdef"), "the test mona lisa.jpg"),
    }
    mock_mime_validation = mocker.patch(
        "sources.services.file_attachments.magic.Magic.from_buffer"
    )
    # we mock the python magic library when checking mimetype
    mock_mime_validation.return_value = "image/jpeg"
    response = client.post(
        "/_upload_experimental_data", data=form_data, content_type="multipart/form-data"
    )
    # assert response is as expected and successful.
    assert response.status_code == 200
    assert len(response.json.get("uploaded_files")) == 1, "1 File should be uploaded"
    assert isinstance(response.json["uploaded_files"][0].get("uuid"), str)
    assert response.json["uploaded_files"][0].get("name") == "the test mona lisa.jpg"
    assert response.json["uploaded_files"][0].get("autogenerated") is False
