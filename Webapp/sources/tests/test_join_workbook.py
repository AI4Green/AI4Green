from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time
from selenium.webdriver.support.select import Select


# test all the buttons are rendered and linked
class JoinWorkBookTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    def apply_to_join_workbook(self, username="SM_Test", password="SM_login"):
        login(self, username, password)
        select_workgroup(self)
        self.driver.find_element_by_id("join-workbook").click()
        select_wb = Select(self.driver.find_element_by_id("workbooks"))
        select_wb.select_by_index(1)
        self.driver.find_element_by_id("submit").click()
        self.assertIn(
            'Your membership has been requested. You will receive a notification when your request has been considered.',
            self.driver.page_source)

    def navigate_to_request_PI(self):
        logout(self)
        login(self, "PI_Test", "PI_login")
        check_notifications(self)
        time.sleep(1)
        self.driver.find_element_by_id('wb1').click()
        time.sleep(1)

    def navigate_to_request_result(self):
        time.sleep(1)
        logout(self)
        login(self, "SM_Test", "SM_login")
        time.sleep(1)
        check_notifications(self)

    def test_join_workbook_deny(self):
        """Tests that a user can ask to join workbook and be denied"""
        self.apply_to_join_workbook()
        self.navigate_to_request_PI()
        self.driver.find_element_by_id('request_deny1').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("This user has not been added to this workbook!", alert.text)
        alert.accept()
        self.navigate_to_request_result()
        self.assertIn("has been denied.", self.driver.page_source)

    def test_join_workbook_approve(self):
        """Tests that a user can ask to join workbook and be denied"""
        self.apply_to_join_workbook()
        self.navigate_to_request_PI()
        self.driver.find_element_by_id('request_approve1').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("This user has been added to this workbook", alert.text)
        alert.accept()
        self.navigate_to_request_result()
        self.assertIn("has been approved.", self.driver.page_source)

    def test_two_users_apply_to_join_workbook(self):
        """
        Test that when 2 users apply to join workbook that the both requests are made and the PI has 2 notifications
        """
        self.apply_to_join_workbook()
        logout(self)
        self.apply_to_join_workbook("SM1_Test", "SM1_login")
        logout(self)
        login(self, "PI_Test", "PI_login")
        self.navigate_to_request_PI()
        self.driver.find_element_by_id("request_approve2").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("This user has been added to this workbook", alert.text)
        alert.accept()
        self.driver.find_element_by_id('request_approve1').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("This user has been added to this workbook", alert.text)
        alert.accept()
        self.navigate_to_request_result()
        self.assertIn("has been approved.", self.driver.page_source)