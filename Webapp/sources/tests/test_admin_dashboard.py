from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db


# test admin dashboard
class AdminDashboardTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    # set up, register and log in
    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        self.login()

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def login(self):
        """This function sets the headless browser up, logs the user in and is called from setUp"""
        setup_selenium(self)
        login(self, "admin", "admin_login")
        self.driver.find_element_by_id('admin_dashboard').click()

    # test admin functions
    def test_admin_functions(self):
        # test admin dashboard has been reached
        self.assertIn("<h2>Admin Dashboard</h2>", self.driver.page_source)
        # check database
        self.assertIn(
            "Users\n                        </td>\n                        <td>\n                            10",
            self.driver.page_source)
        self.assertIn(
            "Workgroups\n                        </td>\n                        <td>\n                            6",
            self.driver.page_source)
        self.assertIn(
            "Compounds\n",
            self.driver.page_source)
        self.assertIn(
            "Reactions\n                        </td>\n                        <td>\n                            2",
            self.driver.page_source)
        # check users
        self.assertIn(
            "PI_Test\n                        </td>\n                        <td>\n                            "
            "PI@test.com\n                        </td>\n                        <td>\n                            "
            "Pat Inglis\n                        </td>\n                        <td>\n                            "
            "Standard",
            self.driver.page_source)
        # check workgroups
        self.assertIn("Test-Workgroup", self.driver.page_source)
