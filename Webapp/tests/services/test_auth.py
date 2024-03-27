from flask.testing import FlaskClient
from tests.utils import login


def test_login(client: FlaskClient):
    """Tests we log in successfully and are redirected to the homepage"""
    login(client)


def test_register(client: FlaskClient):
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
    # now confirm we can log in using these credentials
    login(client, "new_test_user", "new_password")


def test_failed_login(client: FlaskClient):
    """Tests login is unsuccessful with incorrect credentials"""
    login_response = client.post(
        "auth/login",
        data={
            "username": "not_a_registered_user",
            "password": "pastword",
            "submit": True,
        },
    )
    assert (
        login_response.status_code == 302 and login_response.location == "/auth/login"
    ), "login should have failed"

    # Ensure that the user is logged in by checking if the session contains user information
    with client.session_transaction() as session:
        assert "_user_id" not in session, "user should not be logged in"
