from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# test all the buttons are rendered and linked
class JoinWorkgroupTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    # set up, register and log in
    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def apply_to_join_workgroup(self, username="SM_Test", password="SM_login"):
        login(self, username, password)
        self.driver.find_element_by_id("join-WG").click()
        time.sleep(1)
        select_wg = Select(self.driver.find_element_by_id("workgroups"))
        select_wg.select_by_index(1)
        self.driver.find_element_by_id("submit").click()
        self.assertIn(
            'Your membership has been requested. You will receive a notification when your request has been considered.',
            self.driver.page_source)

    def navigate_to_request_PI(self):
        check_notifications(self)
        time.sleep(1)
        self.driver.find_element_by_id('wg1').click()
        time.sleep(1)

    def navigate_to_request_result(self):
        time.sleep(1)
        logout(self)
        login(self, "SM_Test", "SM_login")
        time.sleep(1)
        check_notifications(self)

    def test_join_workgroup_deny(self):
        """Tests that a user can ask to join workgroup and be denied"""
        self.apply_to_join_workgroup()
        logout(self)
        login(self, "PI_Test", "PI_login")
        self.navigate_to_request_PI()
        self.driver.find_element_by_id('request_deny1').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("This request has been denied!", alert.text)
        alert.accept()
        self.navigate_to_request_result()
        self.assertIn("has been denied.", self.driver.page_source)

    def test_join_workgroup_approve(self):
        """Tests that a user can ask to join workgroup and be approved"""
        self.apply_to_join_workgroup()
        logout(self)
        login(self, "PI_Test", "PI_login")
        self.navigate_to_request_PI()
        self.driver.find_element_by_id('request_approve1').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Susan Matthews has been added to Test-Workgroup-2 as a Standard Member!", alert.text)
        alert.accept()
        self.navigate_to_request_result()
        self.assertIn("has been approved!", self.driver.page_source)

    def test_two_users_apply_to_join_workgroup(self):
        """
        Test that when 2 users apply to join workgroup that the both requests are made and the PI has 2 notifications
        """
        self.apply_to_join_workgroup()
        logout(self)
        self.apply_to_join_workgroup("BB_Test", "BB_login")
        logout(self)
        login(self, "PI_Test", "PI_login")
        self.navigate_to_request_PI()
        self.driver.find_element_by_id("request_approve2").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Bob Brown has been added to Test-Workgroup-2 as a Standard Member!", alert.text)
        alert.accept()
        self.driver.find_element_by_id('request_approve1').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Susan Matthews has been added to Test-Workgroup-2 as a Standard Member!", alert.text)
        alert.accept()
        self.navigate_to_request_result()
        self.assertIn("has been approved", self.driver.page_source)

    def test_user_cannot_apply_twice(self):
        """Test when applying a second time that a second notification is not made"""
        self.apply_to_join_workgroup()
        self.driver.find_element_by_id("TopNavHomeButton").click()
        time.sleep(1)
        self.driver.find_element_by_id("join-WG").click()
        time.sleep(1)
        select_wg = Select(self.driver.find_element_by_id("workgroups"))
        select_wg.select_by_index(1)
        self.driver.find_element_by_id("submit").click()

        self.assertIn(
            'You have already submitted a membership request for this workgroup. You will receive a notification '
            'when your request has been considered', self.driver.page_source)


