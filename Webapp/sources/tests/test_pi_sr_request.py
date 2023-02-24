from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# test accepting and denying PI + SR requests
class PIAndSRRequestTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    def request_status(self, username="SR_Test", password="SR_login", requested_status="PI"):
        login(self, username, password)
        select_workgroup(self)
        self.driver.find_element_by_id(f"{requested_status}-status-request").click()
        time.sleep(1)
        self.assertIn(
            'Your request has been made. You will receive a notification when your request has been considered.',
            self.driver.page_source)
        logout(self)

    def handle_request(self, outcome, num=1):
        login(self, "PI_Test", "PI_login")
        check_notifications(self)
        time.sleep(1)
        self.driver.find_element_by_id('wg1').click()
        time.sleep(1)
        self.driver.find_element_by_id(f'request_{outcome}{num}').click()
        time.sleep(1)

    def test_pi_request_deny(self):
        """Tests that a user can ask for pi status and be denied"""
        self.request_status()
        self.handle_request("deny")
        alert = self.driver.switch_to.alert
        self.assertIn("This request has been denied!", alert.text)
        alert.accept()
        time.sleep(1)
        logout(self)
        login(self, "SR_Test", "SR_login")
        time.sleep(1)
        check_notifications(self)
        self.assertIn("has been denied.", self.driver.page_source)

    def test_pi_request_approve(self):
        """Tests that a user can ask for pi status and be denied"""
        self.request_status()
        self.handle_request("approve")
        alert = self.driver.switch_to.alert
        self.assertIn("Sam Reed has changed role from Senior Researcher to Principal Investigator", alert.text)
        alert.accept()
        time.sleep(1)
        logout(self)
        login(self, "SR_Test", "SR_login")
        time.sleep(1)
        check_notifications(self)
        self.assertIn("has been approved", self.driver.page_source)

    def test_sr_request_deny(self):
        """Tests that a user can ask for pi status and be denied"""
        self.request_status("SM_Test", "SM_login", "SR")
        self.handle_request("deny")
        alert = self.driver.switch_to.alert
        self.assertIn("This request has been denied!", alert.text)
        alert.accept()
        time.sleep(1)
        logout(self)
        login(self, "SM_Test", "SM_login")
        time.sleep(1)
        check_notifications(self)
        self.assertIn("has been denied.", self.driver.page_source)

    def test_sr_request_approve(self):
        """Tests that a user can ask for pi status and be denied"""
        self.request_status("SM_Test", "SM_login", "SR")
        self.handle_request("approve")
        alert = self.driver.switch_to.alert
        self.assertIn("Susan Matthews has changed role from Standard Member to Senior Researcher", alert.text)
        alert.accept()
        time.sleep(1)
        logout(self)
        login(self, "SM_Test", "SM_login")
        check_notifications(self)
        self.assertIn("has been approved.", self.driver.page_source)

    def test_2_users_able_to_make_requests(self):
        """Test 2 users can make pi requests"""
        self.request_status()
        self.request_status("SM_Test", "SM_login")

    def test_no_duplicate_requests(self):
        self.request_status()
        login(self, "SR_Test", "SR_login")
        select_workgroup(self)
        self.driver.find_element_by_id(f"PI-status-request").click()
        time.sleep(1)
        self.assertIn(
            'You have already submitted a request for this workgroup. You will receive a notification '
            'when your request has been considered.',
            self.driver.page_source)

    def test_no_user_cannot_request_pi_and_sr_status(self):
        self.request_status("SM_Test", "SM_login")
        login(self, "SM_Test", "SM_login")
        select_workgroup(self)
        self.driver.find_element_by_id(f"SR-status-request").click()
        time.sleep(1)
        self.assertIn(
            'You have already submitted a request for this workgroup. You will receive a notification '
            'when your request has been considered.',
            self.driver.page_source)


