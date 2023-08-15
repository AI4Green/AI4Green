from sources import app
import flask_testing
from database_setup import test_database_create
from auxiliary_for_tests import *


# test the help pages
class ExportDataTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    # set up browser, used to log in
    def setUp(self):
        """This function sets the headless browser up, logs the user in and is called from setUp"""
        setup_selenium(self)
        login(self, "BB_Test", "BB_login")
        select_workgroup(self)

    def test_export_data_button(self):
        """Tests the export data button works"""
        scroll_element(self, "export-csv", True)

    def test_export_pdf(self):
        scroll_element(self, "export-pdf", True)

    def test_export_data_button_before_reactions(self):
        """Tests that the export data button is hidden if there are no reactions"""
        self.driver.find_element_by_id("user-dropdown").click()
        self.driver.find_element_by_id("TopNavLoginButton").click()
        login(self)
        select_workgroup(self)
        self.assertNotIn("""<div id="export-div">""", self.driver.page_source)
