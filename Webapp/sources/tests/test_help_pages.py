from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# test the help pages
class HelpPagesTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    # create the app
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    # set up, register and log in
    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        login(self)
        self.driver.find_element_by_id('InfoButton').click()
        time.sleep(1)

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def test_help_page_video(self):
        """Tests help video button is clickable"""
        self.driver.find_element_by_id("help-video").click()
        self.assertIn('<iframe id="help-video" width="560" height="315"', self.driver.page_source)

    def test_manual_downloads(self):
        """Tests manual is downloadable"""
        # scroll to bottom of page incase footer obscures link
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id("get-manual").click()

    def test_quick_manual_downloads_help_page(self):
        """Tests quick manual is downloadable"""
        # scroll to bottom of page incase footer obscures link
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id("quickstart-guide").click()

    def test_quick_manual_downloads_home_page(self):
        """Tests quick manual is downloadable"""
        self.driver.find_element_by_id("TopNavHomeButton").click()
        time.sleep(1)
        self.driver.find_element_by_id("quickstart-guide").click()

    def test_privacy_link(self):
        """Tests privacy link"""
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id("privacy_link").click()
        # change to the 2nd tab since link opens in new tab
        chwd = self.driver.window_handles[1]
        self.driver.switch_to.window(chwd)
        time.sleep(1)
        self.assertIn('AI4Green Privacy Policy', self.driver.page_source)

    def test_hazard_link(self):
        """Tests privacy link"""
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id("hazard_link").click()
        # change to the 2nd tab since link opens in new tab
        chwd = self.driver.window_handles[1]
        self.driver.switch_to.window(chwd)
        time.sleep(1)
        self.assertIn('Hazard Disclaimer', self.driver.page_source)

    def test_marvin_js_link(self):
        """Tests marvin js link"""
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id("marvin-js-help").click()
        self.assertIn('Marvin JS Help Centre', self.driver.page_source)

    def test_github_link(self):
        """Tests github link"""
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id("github_link").click()
        # change to the 2nd tab since link opens in new tab
        chwd = self.driver.window_handles[1]
        self.driver.switch_to.window(chwd)
        time.sleep(2)
        self.assertIn('AI4Green Installation Guide', self.driver.page_source)
