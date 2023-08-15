from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create
from selenium.common.exceptions import NoAlertPresentException


# test reaction tutorial
class TutorialTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        select_workbook(self)
        make_new_reaction(self)
        self.driver.find_element_by_id('tutorial-mode').click()
        time.sleep(1)
        chwd = self.driver.window_handles[1]
        self.driver.switch_to.window(chwd)
        time.sleep(1)

    def tutorial_proceed(self, tut_id):
        """Click either the tutorial button or alert if present then tutorial button"""
        while True:
            try:
                self.driver.switch_to.alert.dismiss()
            except NoAlertPresentException:
                self.driver.find_element_by_id(tut_id).click()
                time.sleep(1)
                break

    # test if tutorial mode works by clicking through the tutorial
    def test_reaction_tutorial(self):
        # forwards direction
        for i in range(19):
            self.tutorial_proceed("tut-" + str(i + 1) + "-next")
        # check at end
        self.assertIn("Welcome to the reaction sketcher tutorial", self.driver.page_source)
        # backwards direction
        for i in range(19, 0, -1):
            self.tutorial_proceed("tut-" + str(i + 1) + "-back")
        # check back at beginning
        self.assertIn("Congratulations!", self.driver.page_source)
