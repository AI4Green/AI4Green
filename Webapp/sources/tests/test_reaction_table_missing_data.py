from auxiliary_for_tests import *
import flask_testing
from sources import app, db
import time
from database_setup import test_database_create
from selenium import webdriver
from selenium.webdriver.common.keys import Keys


def check_missing_data_number(self, el):
    # check red initially
    element = self.driver.find_element_by_id(el)
    self.assertIn(element.value_of_css_property("border-color"), "rgb(255, 0, 0)")
    # add value
    element.send_keys(5)
    # check now black
    self.assertIn(element.value_of_css_property("border-color"), "rgb(0, 0, 0)")
    # check more data required
    self.driver.find_element_by_id('action-summary').click()
    time.sleep(1)
    alert = self.driver.switch_to.alert
    self.assertIn("Ensure you have entered all the necessary information!", alert.text)
    alert.accept()


def check_missing_data_phys_form(self, el):
    # check red initially
    element = self.driver.find_element_by_id(el)
    self.assertIn(element.value_of_css_property("border-color"), "rgb(255, 0, 0)")
    # try this 3 times - function doesnt always work.
    for i in range(3):
        try:
            actions = webdriver.ActionChains(self.driver)
            actions.move_to_element(element).click()
            time.sleep(1)
            actions.send_keys(Keys.ARROW_DOWN).perform()
            actions.send_keys(Keys.ARROW_DOWN).perform()
            actions.send_keys(Keys.ENTER).perform()
            time.sleep(1)
            # check now black
            self.assertIn(element.value_of_css_property("border-color"), "rgb(0, 0, 0)")
            break
        except:
            continue
    # check more data required
    self.driver.find_element_by_id('action-summary').click()
    time.sleep(1)
    alert = self.driver.switch_to.alert
    self.assertIn("Ensure you have entered all the necessary information!", alert.text)
    alert.accept()


# test missing data
class MissingDataTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    # test missing mass
    def test_missing_data(self):
        demo_reaction(self)
        # check text
        self.assertIn("Please fill in the highlighted boxes to proceed", self.driver.page_source)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        self.driver.find_element_by_class_name('js-add-solvent').click()
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(2)
        # check alert
        alert = self.driver.switch_to.alert
        self.assertIn("Ensure you have entered all the necessary information!", alert.text)
        alert.accept()
        # limiting mass
        check_missing_data_number(self, "js-reactant-rounded-mass1")
        # reactant 1 physical form
        check_missing_data_phys_form(self, "js-reactant-physical-form1")
        # reactant 2 equivalent
        check_missing_data_number(self, "js-reactant-equivalent2")
        # reactant 2 physical form
        check_missing_data_phys_form(self, "js-reactant-physical-form2")
        # reagent 1 equivalent
        check_missing_data_number(self, "js-reagent-equivalent1")
        # reagent 1 physical form
        check_missing_data_phys_form(self, "js-reagent-physical-form1")
        # solvent 1 volume
        check_missing_data_number(self, "js-solvent-volume1")
        # solvent physical form
        check_missing_data_phys_form(self, "js-solvent-physical-form1")
        # product physical form
        check_missing_data_phys_form(self, "js-product-physical-form1")
        # reagent 1 name
        # check red initially
        element = self.driver.find_element_by_id("js-reagent1")
        self.assertIn(element.value_of_css_property("border-color"), "rgb(255, 0, 0)")
        # add value
        actions = webdriver.ActionChains(self.driver)
        element.click()
        time.sleep(1)
        actions.send_keys("ammonia")
        time.sleep(1)
        actions.send_keys(Keys.ENTER).perform()
        time.sleep(1)
        # check now black
        element = self.driver.find_element_by_id("js-reagent1")
        self.assertIn(element.value_of_css_property("border-color"), "rgb(0, 0, 0)")
        # check more data required
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        self.assertIn("Ensure you have entered all the necessary information!", alert.text)
        alert.accept()
        # solvent 1 name
        # check red initially
        element = self.driver.find_element_by_id("js-solvent1")
        self.assertIn(element.value_of_css_property("border-color"), "rgb(255, 0, 0)")
        # add value
        element.click()
        move_down_one_option_action_chains(self)
        # check now black
        # check no more data required
        self.assertIn(element.value_of_css_property("border-color"), "rgb(0, 0, 0)")
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        self.assertIn("Other Risks, controls,", self.driver.page_source)
