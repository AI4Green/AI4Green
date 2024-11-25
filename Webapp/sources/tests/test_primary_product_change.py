from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from selenium.webdriver.support.select import Select
from sources import app, db
import time


# test solvent selection
class PrimaryProductChangeTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        login(self)
        select_workgroup(self)
        select_workbook(self)
        make_new_reaction(self)

    # test products in reaction table
    def test_products_in_reaction_table(self):
        demo_reaction(self)
        pg_source = self.driver.page_source
        self.assertIn('N-Ethylbenzamide', pg_source)
        # self.assertIn('Ethane', pg_source)

    # test the first product is main in reaction and summary tables
    def test_first_product(self):
        demo_reaction(self)
        """First we check that the first product is selected by default in the reaction table"""
        radio_button_selected = self.driver.find_element_by_id('js-main-product1').is_selected()
        self.assertTrue(radio_button_selected)
        """Now we check that the main product has fields for Unreactant and Product Mass"""
        self.driver.find_element_by_id('js-reactant-rounded-mass1').send_keys("500")  # limiting reactant mass
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys(1)
        phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                         'js-product-physical-form1': 1}
        fill_in_physical_forms(self, phys_form_dic)
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        unreacted_field = self.driver.find_element_by_id('js-unreacted-reactant-mass').get_attribute("class")
        product_mass_field = self.driver.find_element_by_id('js-real-product-mass').get_attribute("class")
        self.assertEqual(unreacted_field, "main-product1")
        self.assertEqual(product_mass_field, "main-product1")

    # test the second product is main in reaction and summary tables
    def test_second_product(self):
        add_reaction_sketcher(self, "CC.OC(=O)C1=CC=CC=C1>>CCNC(=O)C1=CC=CC=C1.CC")
        time.sleep(1)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(1)
        """First we select the second product in the reaction table"""
        self.driver.find_element_by_id('js-main-product2').click()
        """Now we check that the main product has fields for Unreactant and Product Mass"""
        self.driver.find_element_by_id('js-reactant-rounded-mass1').send_keys("500")  # limiting reactant mass
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys(1)
        phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                         'js-product-physical-form1': 1, 'js-product-physical-form2': 1}
        fill_in_physical_forms(self, phys_form_dic)
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        unreacted_field = self.driver.find_element_by_id('js-unreacted-reactant-mass').get_attribute("class")
        product_mass_field = self.driver.find_element_by_id('js-real-product-mass').get_attribute("class")
        self.assertEqual(unreacted_field, "main-product2")
        self.assertEqual(product_mass_field, "main-product2")
