from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create



# test all the buttons are rendered and linked
class ReactionsListTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        login(self, "BB_Test", "BB_login")
        select_workgroup(self)

    def test_initially_az(self):
        """Tests that reaction list is initially sorted a to z"""
        self.assertIn('<div id="reaction-1"><h4>a nitro reduction2</h4></div>', self.driver.page_source)
        self.assertIn('<div id="reaction-2"><h4>nitro reduction2</h4></div>', self.driver.page_source)

    def test_change_time(self):
        """Tests that reaction list can be changed to be sorted by time"""
        self.driver.find_element_by_id('time').click()
        time.sleep(1)
        self.assertIn('<div id="reaction-1"><h4>nitro reduction2</h4></div>', self.driver.page_source)
        self.assertIn('<div id="reaction-2"><h4>a nitro reduction2</h4></div>', self.driver.page_source)

    def test_remove_reaction(self):
        """Tests that access to manage workgroup revoked if remove self"""
        self.driver.find_element_by_id("nitro reduction2-delete").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        alert.dismiss()
        time.sleep(1)
        self.assertIn("<h4>nitro reduction2</h4>", self.driver.page_source)
        self.driver.find_element_by_id("nitro reduction2-delete").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        alert.accept()
        time.sleep(1)
        self.assertNotIn("<h4>nitro reduction2</h4>", self.driver.page_source)
