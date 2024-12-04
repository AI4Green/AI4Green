from auxiliary_for_tests import *
from pony.orm import db_session, select
import re
import flask_testing
from database_setup import test_database_create
from sources.email_methods import send_notification_email_test
from sources import app, db


# test notification emails
class NotificationEmailsTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        login(self)

    @db_session
    def test_notification_email(self):
        """Tests notification email shows correct user and the link to notification page works"""
        person = select(x for x in db.Person if x.user.email == "PI@test.com").first()
        email = send_notification_email_test(person)
        # test username is in email
        self.assertIn("Dear PI_Test", email)
        # test link in email
        # find all links in email
        urls = re.findall(r'http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\(\),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+', email)
        for url in urls:
            self.assertIn("http://localhost/notifications", url)
