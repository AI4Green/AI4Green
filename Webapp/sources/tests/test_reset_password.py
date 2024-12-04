from auxiliary_for_tests import *
from pony.orm import select, db_session
import flask_testing
from database_setup import test_database_create
from sources.email_methods import send_password_reset_email_test
from sources import app, db
import time


# test password reset
class ResetPassword(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        self.driver.find_element_by_id("accept_cookies").click()

    def test_reset_password_link(self):
        """Tests the reset password link works"""
        time.sleep(1)
        self.driver.find_element_by_id("reset-password").click()
        self.assertIn("<h2>Reset Password</h2>", self.driver.page_source)

    def test_reset_email_text(self):
        """Tests the reset email text"""
        with db_session:
            user = select(u for u in db.User if "PI@test.com" == u.email).first()
        html_text, token = send_password_reset_email_test(user)
        # email contains username
        self.assertIn("PI_Test", html_text)
        # email contains token
        self.assertIn(token, html_text)

    def test_change_password_page(self):
        """Tests reset password link to change password"""
        # get user
        with db_session:
            user = select(u for u in db.User if "PI@test.com" == u.email).first()
        # get token
        html_text, token = send_password_reset_email_test(user)
        # go to url in email
        self.driver.get("http://localhost:8943/reset_password/" + token)
        # check correct page
        self.assertIn("<h2>Reset Your Password</h2>", self.driver.page_source)
        # put in new password
        password_field = self.driver.find_element_by_id('password')
        password_field.clear()
        password_field.send_keys("newpass")
        # re-enter new password and submit
        password2_field = self.driver.find_element_by_id('password2')
        password2_field.clear()
        password2_field.send_keys("newpass")
        self.driver.find_element_by_id('submit').click()
        time.sleep(3)
        # check message that password has been reset shows and on log in page
        self.assertIn("Your password has been reset", self.driver.page_source)
        self.assertIn("<h2>Sign In</h2>", self.driver.page_source)
        # put in new password
        login(self, "PI_Test", "newpass")
        # check log in is successful
        self.assertIn("<h2>Welcome to AI4Green, PI_Test!</h2>", self.driver.page_source)
