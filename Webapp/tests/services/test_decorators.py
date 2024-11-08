import flask_login
import pytest
from flask.testing import FlaskClient
from flask import Flask, jsonify, url_for
from flask_login import current_user, login_user
from pytest_mock import MockFixture
from typing import List
from tests.utils import login
from sources import services
from sources.decorators import (
    principal_investigator_required,
    principal_investigator_or_senior_researcher_required,
    workbook_member_required,
    workgroup_member_required,
    _is_demo,
    _get_from_request
)


@pytest.mark.parametrize("input_value, search_str, request_args, request_form, expected", [
    (None, "workgroup", {"workgroup": "arg_value"}, {}, "arg_value"),
    (None, "workgroup", {}, {"workgroup": "form_value"}, "form_value"),
    ("default_value", "workgroup", {}, {}, "default_value"),
])
def test_get_from_request(input_value, search_str, request_args, request_form, expected, app):
    """Tests _get_from_request from decorators.py. Uses parameters defined above"""
    with app.test_request_context(query_string=request_args, data=request_form):
        assert _get_from_request(input_value, search_str) == expected


@pytest.mark.parametrize(
    "request_args, request_form, expected", [
        ({"demo": "demo"}, {}, True),
        ({"tutorial": "yes"}, {}, True),
        ({}, {"demo": "demo"}, True),
        ({}, {"tutorial": "yes"}, True),
        ({}, {}, False),
    ]
)
def test_is_demo(request_args, request_form, expected, app):
    """Tests _is_demo from decorators.py. Uses parameters defined above"""
    with app.test_request_context(query_string=request_args, data=request_form):
        assert _is_demo() == expected


def mock_function(workgroup, workbook):
    """Mock function used to test security decorators"""
    return jsonify("Access Granted")


def setup_user_type_mocker(mocker: MockFixture, user_type: str):
    """
    Sets up mocker for sources.services.workgroup.get_user_type function
    Args:
        mocker (MockFixture): Fixture for mocking
        user_type (str): return value for function, "principal_investigator", "senior_researcher" or "standard_member"
    """
    mock_get_user_type = mocker.patch("sources.services.workgroup.get_user_type")
    mock_get_user_type.return_value = user_type


def check_response_success(app: Flask, decorated_function: mock_function, workgroup="Test-Workgroup", workbook=None):
    """
    Checks that the decorator applied to the mock function allows user through
    Args:
        app (Flask): Flask application
        decorated_function (mock_function): instance of mock_function() with applied security decorator
        workgroup (str): name of workgroup
        workbook (str): name of workbook
    """
    with app.test_request_context():
        response = decorated_function(workgroup=workgroup, workbook=workbook)
        assert response.status_code == 200
        assert response.json == "Access Granted"


def check_response_failed(app: Flask, decorated_function: mock_function, workgroup="Test-Workgroup", workbook=None):
    """
    Checks that the decorator applied to the mock function redirects user to homepage
    Args:
        app (Flask): Flask application
        decorated_function (mock_function): instance of mock_function() with applied security decorator
        workgroup (str): name of workgroup
        workbook (str): name of workbook
    """
    with app.test_request_context():
        response = decorated_function(workgroup=workgroup, workbook=workbook)
        assert response.status_code == 302
        assert response.location == url_for("main.index")


def test_principal_investigator_required(client: FlaskClient, mocker: MockFixture, app: Flask):
    """Tests principal_investigator_required decorator"""
    decorated_function = principal_investigator_required(mock_function)

    # test user is PI
    setup_user_type_mocker(mocker, user_type="principal_investigator")
    check_response_success(app, decorated_function)

    # test user is not PI
    setup_user_type_mocker(mocker, user_type="not_PI")
    check_response_failed(app, decorated_function)


def test_principal_investigator_or_senior_researcher_required(client: FlaskClient, mocker: MockFixture, app: Flask):
    """Tests principal_investigator_or_senior_researcher_required decorator"""
    decorated_function = principal_investigator_or_senior_researcher_required(mock_function)

    # mock user is PI
    setup_user_type_mocker(mocker, user_type="principal_investigator")
    check_response_success(app, decorated_function)

    # mock user is SR
    setup_user_type_mocker(mocker, user_type="senior_researcher")
    check_response_success(app, decorated_function)

    # mock user is not PI or SR
    setup_user_type_mocker(mocker, user_type="not_PI")
    check_response_failed(app, decorated_function)


def test_workgroup_member_required(client: FlaskClient, app: Flask):
    """Tests workgroup_member_required decorator"""
    with client:
        login(client)
        decorated_function = workgroup_member_required(mock_function)

        # check default workgroup (user should be a member of the default value)
        check_response_success(app, decorated_function)

        # check a workgroup the user is not a member of
        check_response_failed(app, decorated_function, workgroup="Other-Workgroup")


def test_workbook_member_required(client: FlaskClient,app: Flask):
    """Tests workbook_member_required decorator"""
    with app.app_context():
        with client:
            login(client)
            decorated_function = workbook_member_required(mock_function)

            # check a workbook the user is a member of
            check_response_success(app, decorated_function, workbook="Test-Workbook")

            # check a workbook the user is not a member of
            check_response_failed(app, decorated_function, workbook="Other-Workbook")
