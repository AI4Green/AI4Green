from auxiliary_for_tests import *
from flask_testing import TestCase, LiveServerTestCase
from database_setup import test_database_create
from sources import app, db



# test cookie banner
class CookieBannerTest(LiveServerTestCase, TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        self.start()

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def start(self):
        """This function sets the headless browser up, and takes user to registration page"""
        setup_selenium(self)

    def test_banner_present(self):
        cookie_banner = self.driver.find_element_by_id('cookie-consent-banner').is_displayed()
        self.assertEqual(True, cookie_banner)

    def test_banner_disappears(self):
        self.driver.find_element_by_id('accept_cookies').click()
        cookie_banner = self.driver.find_element_by_id('cookie-consent-banner').is_displayed()
        self.assertEqual(False, cookie_banner)