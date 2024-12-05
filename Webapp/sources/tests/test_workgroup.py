from auxiliary_for_tests import *
from sources import app, db
import flask_testing
from database_setup import test_database_create


# test all the buttons are rendered and linked
class PIButtonsTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        login(self)
        select_workgroup(self)

    def test_workgroup_name_renders(self):
        # test correct workgroup is shown as the title
        self.assertIn("<h2>Test-Workgroup</h2>", self.driver.page_source)

    def test_user_type_renders(self):
        # test correct user type is shown as the title
        self.assertIn("User Type: Principal Investigator", self.driver.page_source)

    def test_change_workgroup_button(self):
        """Tests the change workgroup button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("change-workgroup").click()
        self.assertIn("Welcome to AI4Green, PI_Test!", self.driver.page_source)

    def test_new_reaction_button(self):
        """Tests the new reaction button is clickable and the correct page is loaded"""
        # before a workbook has been check new reaction button not visible
        self.assertIn("""<div id="reaction-content" style="display: block;">""", self.driver.page_source)
        # select workbook and check workbook is transferred to sketcher
        select_workbook(self)
        make_new_reaction(self)
        self.assertIn("Please sketch ", self.driver.page_source)

    def test_manage_workgroup_button(self):
        """Tests the manage workgroup button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("manage-workgroup").click()
        self.assertIn("<h2>Manage Workgroup</h2>", self.driver.page_source)

    def test_manage_workbook_button(self):
        """Tests the manage workbook button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("manage-workbook").click()
        self.assertIn("<h2>Manage Workbooks</h2>", self.driver.page_source)


class SRButtonsTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        login(self, "SR_Test", "SR_login")
        select_workgroup(self)

    def test_workgroup_name_renders(self):
        # test correct workgroup is shown as the title
        self.assertIn("<h2>Test-Workgroup</h2>", self.driver.page_source)

    def test_user_type_renders(self):
        # test correct user type is shown as the title
        self.assertIn("User Type: Senior Researcher", self.driver.page_source)

    def test_change_workgroup_button(self):
        """Tests the change workgroup button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("change-workgroup").click()
        self.assertIn("Welcome to AI4Green, SR_Test!", self.driver.page_source)

    def test_new_reaction_button(self):
        """Tests the new reaction button is clickable and the correct page is loaded"""
        # before a workbook has been check new reaction button not visible
        self.assertIn("""<div id="reaction-content" style="display: block;">""", self.driver.page_source)
        # select workbook and check workbook is transferred to sketcher
        select_workbook(self)
        make_new_reaction(self)
        self.assertIn("Please sketch ", self.driver.page_source)
        self.assertIn("Please sketch ", self.driver.page_source)

    def test_manage_workbook_button(self):
        """Tests the manage workbook button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("manage-workbook").click()
        self.assertIn("<h2>Manage Workbooks</h2>", self.driver.page_source)


class WorkgroupTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        login(self, "SM_Test", "SM_login")
        select_workgroup(self)

    def test_workgroup_name_renders(self):
        # test correct workgroup is shown as the title
        self.assertIn("<h2>Test-Workgroup</h2>", self.driver.page_source)

    def test_user_type_renders(self):
        # test correct user type is shown as the title
        self.assertIn("User Type: Standard", self.driver.page_source)

    def test_change_workgroup_button(self):
        """Tests the change workgroup button is clickable and the correct page is loaded"""
        self.driver.find_element_by_id("change-workgroup").click()
        self.assertIn("Welcome to AI4Green, SM_Test!", self.driver.page_source)

    def test_new_reaction_button(self):
        """Tests the new reaction button is clickable and the correct page is loaded"""
        # before a workbook has been check new reaction button not visible
        self.assertIn("""<div id="reaction-content" style="display: block;">""", self.driver.page_source)
        # select workbook and check workbook is transferred to sketcher
        make_new_reaction(self)
        self.assertIn("Please sketch ", self.driver.page_source)
