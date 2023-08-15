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

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)

    def test_notification_read(self):
        """Tests that a user can ask to join workgroup and be denied"""
        login(self, "SM_Test", "SM_login")
        self.driver.find_element_by_id("join-WG").click()
        time.sleep(1)
        select_wg = Select(self.driver.find_element_by_id("workgroups"))
        select_wg.select_by_index(1)
        self.driver.find_element_by_id("submit").click()
        self.assertIn('Your membership has been requested. You will receive a notification when your request has been considered.', self.driver.page_source)
        logout(self)
        login(self, "PI_Test", "PI_login")
        time.sleep(1)
        # test read notification number shows
        self.assertIn('<span class="badge badge-pill badge-danger">1</span>', self.driver.page_source)
        check_notifications(self)
        time.sleep(1)
        # test notification number no longer shows
        self.assertIn("<sup></sup>", self.driver.page_source)


