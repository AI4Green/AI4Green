import os
from email import generator
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from typing import List

from flask import current_app
from flask_mail import Mail, Message


class EmailSender:
    """
    Custom Flask extension to handle email sending.

    The extension provides capability to switch between sending emails using Flask-Mail, or saving them locally.
    This behaviour depends on the `MAIL_USE_LOCAL` configuration.

    """

    def __init__(self, app=None) -> None:
        """
        Initialise the extension.

        Args:
            app: Flask app instance.
        """
        self.mail = None
        if app is not None:
            self.init_app(app)

    def init_app(self, app) -> None:
        """
        Initialise the Flask-Mail extension with the Flask app.

        Args:
            app: Flask app instance.
        """
        self.mail = Mail(app)

    def send_email(
        self,
        subject: str,
        sender: str,
        recipients: List[str],
        text_body: str,
        html_body: str,
    ) -> None:
        """
        Send an email based on the parameters.

        If config setting MAIL_USE_LOCAL is 'local', the email is saved locally instead.

        Args:
            subject: Subject of the email.
            sender: Sender of the email.
            recipients: Recipients email address.
            text_body: Text body of the email.
            html_body: HTML body of the email.
        """
        if current_app.config["MAIL_USE_LOCAL"] == "local":
            self._save_local(subject, sender, recipients, text_body, html_body)
        else:
            msg = Message(subject, sender=sender, recipients=recipients)
            msg.body = text_body
            msg.html = html_body
            self.mail.send(msg)

    def _save_local(
        self,
        subject: str,
        sender: str,
        recipients: List[str],
        text_body: str,
        html_body: str,
    ) -> None:
        """
        Save the email content to local storage.

        Args:
              subject: Subject of the email.
              sender: Sender of the email.
              recipients: Recipients email address.
              text_body: Text body of the email.
              html_body: HTML body of the email.
        """
        msg = MIMEMultipart()
        msg["To"] = ", ".join(recipients)
        msg["From"] = sender
        msg["Subject"] = subject
        msg.attach(MIMEText(text_body, "plain"))

        # Create the directory if it doesn't exist
        save_dir = current_app.config.get("MAIL_SAVE_DIR")
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        file_path = os.path.join(save_dir, f"{subject} - {recipients[0]}.eml")
        with open(file_path, "w") as f:
            gen = generator.Generator(f)
            gen.flatten(msg)
