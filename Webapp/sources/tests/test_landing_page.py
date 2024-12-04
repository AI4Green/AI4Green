from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time
import unittest


# test all the buttons are rendered and linked
class StandardUserButtonsTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    def test_join_workgroup_button(self):
        """Tests the join workgroup button is clickable and the correct page is loaded"""
        login(self, "SR_Test", "SR_login")
        self.driver.find_element_by_id("join-WG").click()
        time.sleep(1)
        self.assertIn("<h2>Request to Join a Workgroup</h2>", self.driver.page_source)

    def test_create_workgroup_button(self):
        """Tests the create workgroup button is clickable and the correct page is loaded"""
        login(self, "SR_Test", "SR_login")
        self.driver.find_element_by_id("create-workgroup").click()
        time.sleep(1)
        self.assertIn("<h2>Create Workgroup</h2>", self.driver.page_source)

    def test_select_first_workgroup(self):
        """Tests selecting the first workgroup and clicking the proceed to workgroup button"""
        login(self, "SR_Test", "SR_login")
        select_wg = Select(self.driver.find_element_by_id("WG-select"))
        select_wg.select_by_index(1)
        self.driver.find_element_by_id("go-to-workgroup").click()
        self.assertIn("<h2>Test-Workgroup</h2>", self.driver.page_source)

    def test_select_second_workgroup(self):
        """Tests selecting the second workgroup and clicking the proceed to workgroup button"""
        login(self, "SR_Test", "SR_login")
        select_wg = Select(self.driver.find_element_by_id("WG-select"))
        select_wg.select_by_index(2)
        self.driver.find_element_by_id("go-to-workgroup").click()
        self.assertIn("<h2>Test-Workgroup-2</h2>", self.driver.page_source)

    def test_user_info_renders(self):
        login(self, "SR_Test", "SR_login")
        self.assertIn("Welcome to AI4Green, SR_Test!", self.driver.page_source)

    def test_admin_button_hidden(self):
        """Tests non-admin users do not have the admin button"""
        login(self, "SR_Test", "SR_login")
        self.assertNotIn("admin_dashboard", self.driver.page_source)

    def test_admin_link(self):
        login(self, "admin", "admin_login")
        self.assertIn('admin_dashboard', self.driver.page_source)
        self.driver.find_element_by_id("admin_dashboard").click()
        self.assertIn("<h2>Admin Dashboard</h2>", self.driver.page_source)


if __name__ == '__main__':
    unittest.main()
