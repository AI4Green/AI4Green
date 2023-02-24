from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create
from selenium.webdriver.support.select import Select
import unittest


app.config.update(dict(
    TESTING=True,
    DEBUG=False,
    WTF_CSRF_ENABLED=False,
    LOGIN_DISABLED=True
))


# test summary rendered
class SummaryTableRenderTest(flask_testing.LiveServerTestCase, unittest.TestCase):
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

    # test summary table renders
    def test_summary_table_rendered(self):
        """Tests summary table renders"""
        demo_reaction(self)
        time.sleep(1)
        """Let's fill in the reaction table"""
        self.driver.find_element_by_id('js-reactant-rounded-mass1').send_keys("500")
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys("2")
        phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                         'js-product-physical-form1': 1}
        fill_in_physical_forms(self, phys_form_dic)
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        self.assertIn('<b>Reactants/catalysts/reagents</b>', self.driver.page_source)

    def fill_in_demo_mandatory_fields(self):
        time.sleep(2)
        self.driver.find_element_by_id('js-reactant-rounded-mass1').send_keys("500")
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys("2")
        phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                         'js-product-physical-form1': 1}
        fill_in_physical_forms(self, phys_form_dic)

    def test_element_sustainability_green(self):
        """Tests the expected_result is green"""
        demo_reaction(self)
        self.fill_in_demo_mandatory_fields()
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        element_sustainability_color = self.driver.find_element_by_id('js-elements-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(0, 255, 0, 1)', element_sustainability_color)

    def test_element_sustainability_yellow(self):
        """Tests the expected_result is yellow due to Palladium"""
        rxn_smiles = '[Pd].CO>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        self.fill_in_demo_mandatory_fields()
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        element_sustainability_color = self.driver.find_element_by_id('js-elements-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(255, 255, 0, 1)', element_sustainability_color)

    def test_element_sustainability_red(self):
        """Tests the expected_result is red due to Gold"""
        rxn_smiles = 'CO.[Au]>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        self.fill_in_demo_mandatory_fields()
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        element_sustainability_color = self.driver.find_element_by_id('js-elements-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(255, 0, 0, 1)', element_sustainability_color)

