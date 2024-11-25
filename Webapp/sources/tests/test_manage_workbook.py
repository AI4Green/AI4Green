from auxiliary_for_tests import *
from selenium.webdriver.support.select import Select
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# test workbook management
class WorkBookManagementTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        self.driver.find_element_by_id("manage-workbook").click()
        time.sleep(1)

    def test_promote_demote_user(self):
        """Tests that a user can be promoted from user to not user"""
        self.driver.find_element_by_id("remove1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("This user has been removed from the workbook!", alert.text)
        alert.accept()
        self.driver.find_element_by_id("add1").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("This user has been added to this workbook!", alert.text)
        alert.accept()

    def test_change_workbook(self):
        """Tests that a user can change workbook"""
        select_wb = Select(self.driver.find_element_by_id("workbooks"))
        self.assertEqual("Test-Workbook", select_wb.first_selected_option.text)
        select_wb.select_by_index(1)
        select_wb = Select(self.driver.find_element_by_id("workbooks"))
        self.assertEqual("Test-Workbook2", select_wb.first_selected_option.text)
        time.sleep(1)
        self.assertIn('<button id="add1" value="SM@test.com"', self.driver.page_source)
