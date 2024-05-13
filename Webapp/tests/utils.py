from typing import Dict


def login(client, username: str = "test_username", password: str = "test_pw"):
    """Frequently used function to log a user in to gain access to test @login_required functions."""

    login_response = client.post(
        "auth/login", data={"username": username, "password": password, "submit": True}
    )
    assert (
        login_response.status_code == 302 and login_response.location == "/home"
    ), "login failed"

    # Ensure that the user is logged in by checking if the session contains user information
    with client.session_transaction() as session:
        assert "_user_id" in session, "user should be logged in"


def assert_expected_values(actual_data: Dict, expected_data: Dict):
    for key, expected_value in expected_data.items():
        assert (
            actual_data[key] == expected_value
        ), f"{expected_value} not found in data for key: {key}"
