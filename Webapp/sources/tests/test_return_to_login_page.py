from auxiliary_for_tests import *
import flask_testing
from sources import app, db
from database_setup import test_database_create
import unittest

"""The file starts with these functions which help in setting up the unit tests"""


def login(self, username, password):  # login function
    return self.client.post(
        '/auth/login',
        data=dict(username=username, password=password),
        follow_redirects=True
    )


def login_page(self):  # login page
    return self.client.get(
        '/auth/login',
        follow_redirects=True
    )


def registration_page(self):  # registration page
    return self.client.get(
        '/auth/register',
        follow_redirects=True
    )


def reset_password_page(self):  # reset password page
    return self.client.get(
        '/reset_password_request',
        follow_redirects=True
    )


def register(self, username, fullname, email, password, password2):  # registration function
    return self.client.post(
        '/auth/register',
        data=dict(username=username, fullname=fullname, email=email, password=password, password2=password2),
        follow_redirects=True
    )


def logout(self):  # logout function
    return self.client.get(
        '/auth/logout',
        follow_redirects=True
    )


"""The next tests use the backend to check that a user is redirected to login page"""


class TestReturnToLogin(flask_testing.TestCase):
    """These functions test the response from the login page."""

    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        login(self, 'PI_Test', 'PI_login')

    def tearDown(self):
        logout(self)
        restore_db()

    def test_logout_to_login(self):
        """Test that we return to the login page after logging out"""
        response = logout(self)
        self.assertIn(b'Please log in to access this page', response.data)

    def test_registration_to_login(self):
        """Test that we return to the login page from the registration page"""
        logout(self)
        """First, check that we are on the registration page"""
        response_registration_page = registration_page(self)
        self.assertIn(b'Repeat Password', response_registration_page.data)
        """Second, check that we are back to the login page"""
        response_login_page = login_page(self)
        self.assertIn(b'Sign In', response_login_page.data)

    def test_reset_password_to_login(self):
        """Test that we return to the login page from the reset password page"""
        logout(self)
        """First, check that we are on the reset password page"""
        response_reset_password_page = reset_password_page(self)
        self.assertIn(b'Reset Password', response_reset_password_page.data)
        """Second, check that we are back to the login page"""
        response_login_page = login_page(self)
        self.assertIn(b'Sign In', response_login_page.data)


# class TestBrowserReturnsToLogin(flask_testing.LiveServerTestCase, unittest.TestCase):
#
#     def create_app(self):
#         # Default port is 5000
#         app.config['LIVESERVER_PORT'] = 8943
#         # Default timeout is 5 seconds
#         app.config['LIVESERVER_TIMEOUT'] = 10
#         db.drop_all_tables(with_all_data=True)
#         db.create_tables()
#         return app
#
#     def browser_login(self):
#         """This function sets the headless browser up, logs the user in and is called from setUp"""
#         op = webdriver.ChromeOptions()
#         op.add_argument('headless')
#         chrome_dir = os.path.join(basedir, "sources/tests/test_data/chromedriver")
#         self.driver = webdriver.Chrome(chrome_dir, options=op)
#         self.driver.get("http://localhost:8943")
#
#         username_field = self.driver.find_element_by_id('username')
#         username_field.clear()
#         username_field.send_keys("tom")
#
#         password_field = self.driver.find_element_by_id('password')
#         password_field.clear()
#         password_field.send_keys("sun")
#
#         self.driver.find_element_by_id('submit').click()
#
#     def setUp(self):
#         """We load a test instance of the app, clear the db, register a new user, then login"""
#         self.app = app.test_client()
#         db.drop_all_tables(with_all_data=True)
#         db.create_tables()
#         register(self, 'tom', 'Tom Waits', 'tom@mail.com', 'sun', 'sun')
#         self.browser_login()
#
#     def test_browser_logout_to_login(self):
#         """This function tests that the browser logs out"""
#         self.driver.find_element_by_link_text('Logout').click()
#         pg_source = self.driver.page_source
#         self.assertIn('Please log in to access this page', pg_source, "The browser has not logged out")
#
#     # def test_marvin_demo_function(self):
#     #     """Here we test that compounds can be loaded into marvin_js by the example demo button by pressing the submit
#     #      button and checking the correct compound appears in the reaction table"""
#     #     pg_source = self.driver.page_source
#     #     self.assertNotIn('benzoic acid', pg_source, "no compound names should be present at this stage")
#     #     time.sleep(3)
#     #     self.driver.find_element_by_id('demo-button').click()
#     #     time.sleep(1)
#     #     self.driver.find_element_by_id('action-button-submit').click()
#     #     time.sleep(15)
#     #     pg_source = self.driver.page_source
#     #     self.assertIn('benzoic acid', pg_source, "benzoic acid should now be loaded in the reaction table")
#
#     def tearDown(self):
#         # close the browser window
#         self.driver.quit()
#
#     @classmethod
#     def tearDownClass(cls):
#         db.drop_all_tables(with_all_data=True)
#         db.create_tables()


if __name__ == '__main__':
    unittest.main()
