from auxiliary_for_tests import *
from sources import app, db
import time
from selenium.webdriver.support.select import Select
import flask_testing
from database_setup import test_database_create


class UnknownReactantTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """ Need to reload the reaction in the sketcher using the data from the reaction database"""

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
        select_workgroup(self)
        select_workbook(self)
        make_new_reaction(self)

    def test_unknown_reactant(self):
        """ This function returns the smiles string of the reaction currently in the chemical sketcher """
        # switch to the marvin js frame
        self.driver.switch_to.frame(0)
        # click on the import button (2nd along) to open import structure popup
        elem_ls = self.driver.find_elements_by_xpath('//*[local-name() = "svg"]')
        elem_ls[1].click()
        # change the format to smiles
        format_dropdown = Select(self.driver.find_elements_by_xpath('//*[@class="gwt-ListBox"]')[0])
        format_dropdown.select_by_visible_text('SMILES')
        # find the textarea in the popup window and save the smiles string as a variable
        text_area = self.driver.find_elements_by_xpath('//*[@class="gwt-TextArea"]')[0]
        text_area.send_keys("C1NCNCN1>>C1CCCCC1")
        # buttons = self.driver.find_element_by_class_name("gwt-Button.mjs-marginTop")
        buttons = self.driver.find_element_by_xpath('//*[@class="gwt-Button mjs-marginTop"]')
        buttons.click()
        time.sleep(1)
        self.driver.switch_to_default_content()
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Reactant 1 not in database', pg_source)

    def test_unknown_product(self):
        """ This function returns the smiles string of the reaction currently in the chemical sketcher """
        # switch to the marvin js frame
        self.driver.switch_to.frame(0)
        # click on the import button (2nd along) to open import structure popup
        elem_ls = self.driver.find_elements_by_xpath('//*[local-name() = "svg"]')
        elem_ls[1].click()
        # change the format to smiles
        format_dropdown = Select(self.driver.find_elements_by_xpath('//*[@class="gwt-ListBox"]')[0])
        format_dropdown.select_by_visible_text('SMILES')
        # find the textarea in the popup window and save the smiles string as a variable
        text_area = self.driver.find_elements_by_xpath('//*[@class="gwt-TextArea"]')[0]
        text_area.send_keys("C1CCCCC1>>C1NCNCN1")
        # buttons = self.driver.find_element_by_class_name("gwt-Button.mjs-marginTop")
        buttons = self.driver.find_element_by_xpath('//*[@class="gwt-Button mjs-marginTop"]')
        buttons.click()
        time.sleep(1)
        self.driver.switch_to_default_content()
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Product 1 not in database', pg_source)
