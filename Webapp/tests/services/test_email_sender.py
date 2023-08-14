import pytest
from flask import current_app
from sources.extensions import mail


def test_send_email_sending(client):
    """
    Test sending an email when MAIL_USE_LOCAL is set to False.

    This test ensures that the 'send_email' function correctly sends an email when the 'MAIL_USE_LOCAL' configuration
    is set to False, indicating that the email should be sent using the Flask-Mail extension.

    Args:
        client: Flask test client.

    """
    current_app.config["MAIL_USE_LOCAL"] = False
    subject = "Test Email"
    sender = "sender@example.com"
    recipients = ["recipient@example.com"]
    text_body = "Text body of the email."
    html_body = "<p>HTML body of the email.</p>"

    with mail.mail.record_messages() as outbox:
        mail.send_email(subject, sender, recipients, text_body, html_body)
        assert len(outbox) == 1
        assert outbox[0].subject == subject


def test_send_email_local_save(client):
    """
    Test saving an email locally when MAIL_USE_LOCAL is set to True.

    This test ensures that the 'send_email' function correctly saves the email locally when the 'MAIL_USE_LOCAL'
    configuration is set to True. The email content should be stored in a file on the local storage.

    Args:
        client: Flask test client.

    """
    subject = "Test Email"
    sender = "sender@example.com"
    recipients = ["recipient@example.com"]
    text_body = "Text body of the email."
    html_body = "<p>HTML body of the email.</p>"

    mail.send_email(subject, sender, recipients, text_body, html_body)

    with open(f"temp/{subject} - {recipients[0]}.eml", "r") as file:
        content = file.read()
        assert subject in content
        assert sender in content
        assert recipients[0] in content
        assert text_body in content
