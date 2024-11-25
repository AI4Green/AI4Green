from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create
import unittest


# test all the buttons are rendered and linked
class StandardUserButtonsTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    # set up, register and log in
    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        self.browser_login()

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    # set up browser, used to log in
    def browser_login(self):
        """This function sets the headless browser up, logs the user in and is called from setUp"""
        setup_selenium(self)
        login(self, "SR_Test", "SR_login")

    # test the Home button is linked
    def test_home_button_linked(self):
        self.driver.find_element_by_id('TopNavHomeButton').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Please select a Workgroup to get started', pg_source)  # check that we are on the Home page

    # test the Logout button is linked
    def test_logout_button_linked(self):
        self.driver.find_element_by_id('user-dropdown').click()
        self.driver.find_element_by_id('TopNavLoginButton').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Sign In', pg_source)  # check that we are on the Login page

    # test the Build Reaction button is linked
    def test_demo_button_linked(self):
        self.driver.find_element_by_id('TopNavSketcherButton').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Please sketch or upload your', pg_source)  # check that we are on the Build Reaction page

    # test the manage account button is linked
    def test_manage_account_button_linked(self):
        self.driver.find_element_by_id('user-dropdown').click()
        self.driver.find_element_by_id('manage-account').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Manage Account', pg_source)  # check that we are on the Update username page


if __name__ == '__main__':
    unittest.main()
