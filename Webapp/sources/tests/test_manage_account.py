from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time

# test account management
class AccountManagementTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        self.driver.find_element_by_id('user-dropdown').click()
        self.driver.find_element_by_id('manage-account').click()
        time.sleep(1)

    def test_change_password_button(self):
        """Tests the change password button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("update-email-password").click()
        self.assertIn("<h2>Update Email or Password</h2>", self.driver.page_source)

    def test_delete_profile_button(self):
        """Tests the delete profile button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("delete-profile").click()
        self.assertIn("Delete Profile", self.driver.page_source)

    def test_delete_profile_confirmation(self):
        """Tests the confirm delete profile button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("delete-profile").click()
        self.driver.find_element_by_id("confirm-delete-profile").click()
        self.assertIn("<h2>Delete Profile Confirmation</h2>", self.driver.page_source)

    def test_delete_profile_pi_workgroup(self):
        """Tests that workgroups which would have no PI are flagged and button works"""
        self.driver.find_element_by_id("delete-profile").click()
        self.driver.find_element_by_id("confirm-delete-profile").click()
        # flag orphaned workgroups
        self.assertIn("Assign New PI", self.driver.page_source)
        # button to workgroup management works
        self.driver.find_element_by_id("assign-new-pi1").click()
        self.assertIn("Manage Workgroup", self.driver.page_source)
        self.assertIn("Test-Workgroup</b>", self.driver.page_source)

    def test_delete_profile_deletion(self):
        """Tests that workgroups which would have no PI are flagged and button works"""
        self.driver.find_element_by_id("delete-profile").click()
        self.driver.find_element_by_id("confirm-delete-profile").click()
        # test profile deletes
        self.driver.find_element_by_id("confirm-delete-profile-delete").click()
        time.sleep(1)
        self.assertIn("Please log in to access this page.", self.driver.page_source)
        self.assertIn("<h2>Sign In</h2>", self.driver.page_source)
        login(self)
        self.assertIn("<h2>Sign In</h2>", self.driver.page_source)

    def test_user_info_renders(self):
        self.assertIn("Username: PI_Test", self.driver.page_source)
        self.assertIn("Email: PI@test.com", self.driver.page_source)
