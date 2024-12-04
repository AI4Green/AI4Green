from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# test reaction locks
class CompleteReactionLockedTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        select_workbook(self, idx=0)
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        time.sleep(8)
        clear_and_send_keys(self, 'js-unreacted-reactant-mass', '1')
        clear_and_send_keys(self, 'js-real-product-mass', '9')
        self.driver.find_element_by_id('complete-reaction-button').click()
        time.sleep(1)

    # test reaction locks immediately
    def test_reaction_locks_immediately(self):
        rxn_save_indicator_text = self.driver.find_element_by_id('reaction-saved-indicator').get_attribute('innerHTML')
        self.assertEqual("Reaction Changes Saved &amp; Locked", rxn_save_indicator_text)
        self.assertIn('''js-reactant-rounded-masses remove-highlight-filled-cell" disabled''', self.driver.page_source)
        self.assertIn('''"js-supervisor" type="text" size="25em" disabled''', self.driver.page_source)

    # test reloaded reaction is still locked
    def test_reaction_reloaded_locked(self):
        self.driver.find_element_by_id('TopNavHomeButton').click()
        select_workgroup(self)
        select_workbook(self, idx=0)
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        time.sleep(8)
        self.assertIn('''js-reactant-rounded-masses remove-highlight-filled-cell" disabled''', self.driver.page_source)
        self.assertIn('''"js-supervisor" type="text" size="25em" disabled''', self.driver.page_source)
