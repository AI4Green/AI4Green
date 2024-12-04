from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
import time
from sources import app, db


# test add new reagent by cas
class AddNewReagentCASTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    # set up, register and log in
    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        login(self, "SR_Test", "SR_login")
        select_workgroup(self)
        select_workbook(self)
        make_new_reaction(self)
        demo_reaction(self)

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    # test if cas not in database add new reagent is opened
    def test_reagent_cas_not_in_database(self):
        self.assertEqual(self.driver.find_element_by_id('js-input-reagent1').value_of_css_property("display"), "none")
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_id('js-reagent1').send_keys('17912-87-7')
        time.sleep(3)
        # alert = self.driver.switch_to.alert
        # self.assertIn("CAS 17912-87-7 is not found in our database", alert.text)
        # alert.accept()
        # time.sleep(1)
        self.assertEqual(self.driver.find_element_by_id('js-input-reagent1').value_of_css_property("display"), "table-row")
