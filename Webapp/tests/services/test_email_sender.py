import os

import pytest
from flask import current_app
from sources.extensions import mail

IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


def test_send_email_sending(app):
    """
    Test sending an email when MAIL_USE_LOCAL is set to False.

    This test ensures that the 'send_email' function correctly sends an email when the 'MAIL_USE_LOCAL' configuration
    is set to False, indicating that the email should be sent using the Flask-Mail extension.

    Args:
        client: Flask test client.

    """
    subject = "Test Email"
    sender = "sender@example.com"
    recipients = ["recipient@example.com"]
    text_body = "Text body of the email."
    html_body = "<p>HTML body of the email.</p>"

    with app.app_context():
        current_app.config["MAIL_USE_LOCAL"] = False
        with mail.mail.record_messages() as outbox:
            mail.send_email(subject, sender, recipients, text_body, html_body)
            assert len(outbox) == 1
            assert outbox[0].subject == subject


@pytest.mark.skipif(IN_GITHUB_ACTIONS, reason="Local temp files don't work in CI")
def test_send_email_local_save(app):
    """
    Test saving an email locally when MAIL_USE_LOCAL is set to 'local'.

    This test ensures that the 'send_email' function correctly saves the email locally when the 'MAIL_USE_LOCAL'
    configuration is set to 'local'. The email content should be stored in a file on the local storage.

    Args:
        client: Flask test client.

    """
    subject = "Test Email"
    sender = "sender@example.com"
    recipients = ["recipient@example.com"]
    text_body = "Text body of the email."
    html_body = "<p>HTML body of the email.</p>"

    with app.app_context():
        mail.send_email(subject, sender, recipients, text_body, html_body)

        file_path = os.path.join(
            app.config["MAIL_SAVE_DIR"], f"{subject} - {recipients[0]}.eml"
        )
        with open(file_path, "r") as file:
            content = file.read()
            assert subject in content
            assert sender in content
            assert recipients[0] in content
            assert text_body in content
