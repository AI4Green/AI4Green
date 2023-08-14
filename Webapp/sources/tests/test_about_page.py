from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# test the help pages
class AboutPageTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        time.sleep(1)

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def test_learn_more_home_button(self):
        """Tests the learn more button works"""
        self.driver.find_element_by_id("about-button").click()
        self.assertIn('<h1>About AI4Green</h1>', self.driver.page_source)

    def test_github_link(self):
        """Tests github link"""
        self.driver.find_element_by_id("about-button").click()
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id("github_link").click()
        # change to the 2nd tab since link opens in new tab
        chwd = self.driver.window_handles[1]
        self.driver.switch_to.window(chwd)
        time.sleep(2)
        self.assertIn('AI4Green Installation Guide', self.driver.page_source)

    def test_help_link(self):
        """Tests help link"""
        self.driver.find_element_by_id("about-button").click()
        # scroll to bottom of page incase footer obscures link
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id('about-to-help').click()
        self.assertIn('AI4Green Documentation', self.driver.page_source)

    def test_email_link(self):
        """Tests email link"""
        self.driver.find_element_by_id("about-button").click()
        # scroll to bottom of page incase footer obscures link
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.assertIn('admin@ai4green.app', self.driver.page_source)
        self.driver.find_element_by_id("email-link").click()

    def test_paper_link(self):
        """Tests paper link"""
        self.driver.find_element_by_id("about-button").click()
        # scroll to bottom of page incase footer obscures link
        self.driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        self.driver.find_element_by_id('paper-link').click()
        # change to the 2nd tab since link opens in new tab
        chwd = self.driver.window_handles[1]
        self.driver.switch_to.window(chwd)
        time.sleep(2)
        self.assertIn('AI4Green: An Open-Source ELN', self.driver.page_source)
