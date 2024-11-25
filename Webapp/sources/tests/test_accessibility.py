from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app
import time



# test colour changing accessibility features
class AccessibilityTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        self.driver.find_element_by_id('user-dropdown').click()
        self.driver.find_element_by_id('accessibility').click()
        time.sleep(1)

    def make_and_fill_in_reaction(self):
        select_workgroup(self)
        select_workbook(self)
        make_new_reaction(self)
        demo_reaction(self)
        """Let's fill in the reaction table"""
        self.driver.find_element_by_id('js-reactant-rounded-mass1').send_keys("500")
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys("2")
        phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                         'js-product-physical-form1': 1}
        fill_in_physical_forms(self, phys_form_dic)
        # check solvent highly hazardous colour correct
        # get benzene in solvent box
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        self.driver.find_element_by_id('js-solvent1').click()
        self.driver.find_element_by_id('js-solvent1').send_keys("ben")
        time.sleep(1)
        move_down_one_option_action_chains(self)

    def test_changing_colours(self):
        """Tests that changing the colours changes the colours in the reaction table, summary table, and solvent flash
        cards"""
        # change the colours
        self.driver.execute_script("$('#acceptable-colours').val('#008000')")
        self.driver.execute_script("$('#acceptable-colours-text').val('#ffffff')")
        self.driver.execute_script("$('#warning-colours').val('#d2b48c')")
        self.driver.execute_script("$('#warning-colours-text').val('#ffffff')")
        self.driver.execute_script("$('#hazard-colours').val('#ff0000')")
        self.driver.execute_script("$('#hazard-colours-text').val('#ffffff')")
        self.driver.execute_script("$('#highly-hazard-colours').val('#cb4154')")
        self.driver.execute_script("$('#highly-hazard-colours-text').val('#ffffff')")
        self.driver.find_element_by_id("change-colours").click()
        time.sleep(1)
        self.driver.switch_to.alert.accept()
        self.driver.find_element_by_id("TopNavHomeButton").click()
        self.make_and_fill_in_reaction()
        solvent_color = self.driver.find_element_by_id('js-solvent1').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(203, 65, 84, 1)', solvent_color)
        solvent_color = self.driver.find_element_by_id('js-solvent1').value_of_css_property(
            'color')
        self.assertEqual('rgba(255, 255, 255, 1)', solvent_color)
        self.driver.find_element_by_class_name('js-remove-solvent').click()
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(2)
        # check summary table elements correct colour
        element_sustainability_color = self.driver.find_element_by_id('js-elements-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(0, 128, 0, 1)', element_sustainability_color)
        element_sustainability_color = self.driver.find_element_by_id('js-elements-cell').value_of_css_property(
            'color')
        self.assertEqual('rgba(255, 255, 255, 1)', element_sustainability_color)
        element_safety_color = self.driver.find_element_by_id('js-safety-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(255, 0, 0, 1)', element_safety_color)
        element_safety_color = self.driver.find_element_by_id('js-safety-cell').value_of_css_property(
            'color')
        self.assertEqual('rgba(255, 255, 255, 1)', element_safety_color)
        element_ae_color = self.driver.find_element_by_id('js-ae-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(210, 180, 140, 1)', element_ae_color)
        element_ae_color = self.driver.find_element_by_id('js-ae-cell').value_of_css_property(
            'color')
        self.assertEqual('rgba(255, 255, 255, 1)', element_ae_color)

    def test_reset_colours(self):
        """Tests that resetting the colors works"""
        self.driver.find_element_by_id("reset-colours").click()
        time.sleep(1)
        self.driver.switch_to.alert.accept()
        self.driver.find_element_by_id("TopNavHomeButton").click()
        self.make_and_fill_in_reaction()
        solvent_color = self.driver.find_element_by_id('js-solvent1').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(139, 0, 0, 1)', solvent_color)
        solvent_color = self.driver.find_element_by_id('js-solvent1').value_of_css_property(
            'color')
        self.assertEqual('rgba(255, 255, 255, 1)', solvent_color)
        self.driver.find_element_by_class_name('js-remove-solvent').click()
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(2)
        # check summary table elements correct colour
        element_sustainability_color = self.driver.find_element_by_id('js-elements-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(0, 255, 0, 1)', element_sustainability_color)
        element_sustainability_color = self.driver.find_element_by_id('js-elements-cell').value_of_css_property(
            'color')
        self.assertEqual('rgba(0, 0, 0, 1)', element_sustainability_color)
        element_safety_color = self.driver.find_element_by_id('js-safety-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(255, 0, 0, 1)', element_safety_color)
        element_safety_color = self.driver.find_element_by_id('js-safety-cell').value_of_css_property(
            'color')
        self.assertEqual('rgba(0, 0, 0, 1)', element_safety_color)
        element_ae_color = self.driver.find_element_by_id('js-ae-cell').value_of_css_property(
            'background-color')
        self.assertEqual('rgba(255, 255, 0, 1)', element_ae_color)
        element_ae_color = self.driver.find_element_by_id('js-ae-cell').value_of_css_property(
            'color')
        self.assertEqual('rgba(0, 0, 0, 1)', element_ae_color)
