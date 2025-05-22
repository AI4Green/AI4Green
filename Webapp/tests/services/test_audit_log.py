from flask import Flask
from flask.testing import FlaskClient
from sources import models
from sources.extensions import db
from sqlalchemy import update
from tests.utils import login


def test_data_access_history(app: Flask, client: FlaskClient):
    """Tests get_data_access_history can retrieve from DataAccessHistory table"""
    login(client)

    # Send a GET request to the route with the form data and confirm response
    response = client.get("/data_access_history")
    # Assert that the response is as expected
    assert response.status_code == 200
    assert response.json[0]["workgroup"] == "Test-Workgroup"
    assert response.json[0]["old_role"] == "Senior Researcher"


def test_data_export_history(app: Flask, client: FlaskClient):
    """Tests get_data_export_history can retrieve from DataExportRequest tables"""
    login(client)

    with app.app_context():
        create_data_export()

    # Send a GET request to the route with the form data and confirm response
    response = client.get("/data_export_history")

    # Assert that the response is as expected
    assert response.status_code == 200
    assert response.json[0]["workgroup"] == "Test-Workgroup"
    assert response.json[0]["workbook"] == "Test-Workbook"
    assert response.json[0]["reactions"] == ["TW1-001"]


def create_data_export():
    """updates the data export request to be approved"""
    update_query = (
        update(models.data_export_request_approvers)
        .where(models.data_export_request_approvers.c.data_export_request_id == 1)
        .where(models.data_export_request_approvers.c.person_id == 1)
        .values(approved=True, responded=True)
    )
    db.session.execute(update_query)
    db.session.commit()
