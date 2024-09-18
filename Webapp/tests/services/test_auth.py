from bokeh.models import Markup
from flask.testing import FlaskClient
from flask import Flask
from tests.utils import login, login_response
from sources import services, db

def test_login(client: FlaskClient):
    """Tests we log in successfully and are redirected to the homepage"""
    login(client)


def test_register(client: FlaskClient, app: Flask):
    """Tests we can register a new user successfully"""
    response = client.post(
        "auth/register",
        data={
            "username": "new_test_user",
            "email": "new_test_user@mail.com",
            "fullname": "George Test",
            "password": "new_password",
            "password2": "new_password",
            "privacy": True,
            "hazard_disclaimer": True,
            "submit": True,
        },
    )

    # upon success user is redirected to login page
    assert response.status_code == 302 and response.location == "/auth/login"

    # users should be notified of email validation
    with client.session_transaction() as session:
        assert  "An email has been sent to your address. Please follow the instructions to verify your account." in session["_flashes"][0], "email validation was not requested"

    # mock email verification to test login with new credentials
    with app.app_context():
        user = services.user.from_email("new_test_user@mail.com")
        user.is_verified = True
        db.session.commit()

    # test login with new credentials
    login(client, username="new_test_user", password="new_password")



def test_failed_login_invalid_credentials(client: FlaskClient):
    """Tests login is unsuccessful with incorrect credentials"""

    response = login_response(client, username="not_a_registered_user", password="pastword")

    # upon success user is redirected to login page
    assert response.status_code == 302 and response.location == "/home", "user should be redirected"

    # Ensure that the user is logged in by checking if the session contains user information
    with client.session_transaction() as session:
        assert "_user_id" not in session, "user should not be logged in"


def test_failed_login_not_validated(client: FlaskClient):
    """Tests login is unsuccessful if user has not validated email"""
    response = login_response(client, username="not_verified", password="not_verified")

    assert response.status_code == 302 and response.location == "/home", "user should be redirected"

    with (client.session_transaction() as session):
        assert "_user_id" not in session, "user should not be logged in"
        assert "Please verify your email address before logging in. Didn't receive an email? <a href='/email_verification_request/2' class='alert-link'>Click here</a> to resend." in session["_flashes"][0], "user was not notified "
