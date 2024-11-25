from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create
from selenium.webdriver.common.keys import Keys


# test reagent selection
class ReagentSelectionTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    # test reagent data returned
    def test_reagent_list(self):
        demo_reaction(self)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_id('js-reagent1').click()
        clear_and_send_keys(self, 'js-reagent1', 'urani')
        actions = webdriver.ActionChains(self.driver)
        actions.send_keys(Keys.ENTER)
        time.sleep(2)
        reagent_option = self.driver.find_element_by_css_selector("option[value='Uranium']").get_attribute("value")
        self.driver.find_element_by_id('js-reagent1').click()
        time.sleep(1)
        self.assertEqual("Uranium", reagent_option)
        reagent_option = self.driver.find_element_by_css_selector("option[value='Uranium dioxide']").\
            get_attribute("value")
        time.sleep(1)
        self.assertEqual("Uranium dioxide", reagent_option)
        reagent_option = self.driver.find_element_by_css_selector("option[value='Uranium nitrate oxide (UO2(NO3)2)']").\
            get_attribute("value")
        time.sleep(1)
        self.assertEqual("Uranium nitrate oxide (UO2(NO3)2)", reagent_option)
