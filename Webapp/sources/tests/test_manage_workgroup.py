from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# test workgroup management
class WorkGroupManagementTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        select_workgroup(self)
        self.driver.find_element_by_id("manage-workgroup").click()
        time.sleep(1)

    def test_promote_demote_user(self):
        """Tests that a user can be promoted from sm to sr and sr to pi"""
        self.driver.find_element_by_id("sr-to-pi1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Sam Reed has changed role from Senior Researcher to Principal Investigator", alert.text)
        alert.accept()
        self.driver.find_element_by_id("sm-to-sr1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Susan Matthews has changed role from Standard Member to Senior Researcher", alert.text)
        alert.accept()
        self.driver.find_element_by_id("pi-to-sr2").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Sam Reed has changed role from Principal Investigator to Senior Researcher", alert.text)
        alert.accept()
        self.driver.find_element_by_id("sr-to-sm1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Susan Matthews has changed role from Senior Researcher to Standard Member", alert.text)
        alert.accept()

    def test_remove_user(self):
        """Tests that a user can be removed"""
        self.driver.find_element_by_id("sm-remove1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Are you sure you want to remove this user entirely", alert.text)
        alert.accept()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Stewart Mathers has been removed from the Test-Workgroup and associated Workbooks", alert.text)
        alert.accept()

    def test_always_pi(self):
        """Tests that sole pi cannot be removed"""
        self.driver.find_element_by_id("pi-remove1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Are you sure you want to remove this user entirely", alert.text)
        alert.accept()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Cannot remove PI", alert.text)
        alert.accept()

    def test_remove_self(self):
        """Tests that access to manage workgroup revoked if remove self"""
        self.driver.find_element_by_id("sr-to-pi1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        alert.accept()
        time.sleep(1)
        self.driver.find_element_by_id("pi-remove1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        alert.accept()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        alert.accept()
        time.sleep(1)
        self.assertIn('You do not have permission to view this page', self.driver.page_source)
