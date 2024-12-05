from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create


# test the help pages
class WBReactionsTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        time.sleep(1)
        select_workgroup(self)

    def test_reactions_from_other_users_visible(self):
        """Tests that reactions created in the same workbook by different users are viewable by each other"""
        make_new_reaction(self, "PI reaction")
        demo_reaction(self)
        mass_field = self.driver.find_element_by_id('js-reactant-rounded-mass1')
        mass_field.clear()
        mass_field.send_keys("500")
        time.sleep(1)
        self.driver.find_element_by_id('js-reactant-equivalent2').click()
        time.sleep(2)
        logout(self)
        login(self, "SR_Test", "SR_login")
        time.sleep(1)
        select_workgroup(self)
        # test it says locked
        self.assertIn('Read-only reaction created by PI_Test', self.driver.page_source)
        # test no delete button
        self.assertNotIn('PI reaction-delete', self.driver.page_source)
        # test reloads but locked
        self.driver.find_element_by_id("PI reaction-reload").click()
        time.sleep(2)
        element = self.driver.find_element_by_id("js-reaction-name2")
        self.assertEqual(element.is_enabled(), False)


