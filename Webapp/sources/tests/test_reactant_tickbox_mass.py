from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create



class LimitingReactantTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        setup_selenium(self)
        login(self, "SR_Test", "SR_login")
        select_workgroup(self)
        select_workbook(self)
        make_new_reaction(self)

    def tearDown(self):
        restore_db()
        self.driver.quit()

    # test reaction table displays reactants
    def test_reactants_in_reaction_table(self):
        demo_reaction(self)
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Ethylamine', pg_source)
        self.assertIn('Benzoic acid', pg_source)

    # test the first reactant is limiting by default
    def test_first_reactant_limiting(self):
        demo_reaction(self)
        time.sleep(1)
        radio_button_selected = self.driver.find_element_by_id('js-reactant-limiting1').is_selected()
        self.assertTrue(radio_button_selected)

    # test the limiting reactant mass is changed following radiobutton
    def test_limiting_reactant_mass(self):
        demo_reaction(self)
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys("2")
        self.driver.find_element_by_id('js-reactant-limiting2').click()  # we set the second reactant as limiting
        """Now we check that the first reactant mass is readonly and the second reactant mass is not"""
        first_reactant_mass_field = self.driver.\
            find_element_by_id('js-reactant-rounded-mass1').get_attribute("readonly")
        second_reactant_mass_field = self.driver.\
            find_element_by_id('js-reactant-rounded-mass2').get_attribute("readonly")
        self.assertTrue(first_reactant_mass_field)
        self.assertFalse(second_reactant_mass_field)
        """Now we check that the second reactant equivalent is readonly and the first reactant equivalent is not"""
        first_reactant_equivalent_field = self.driver.\
            find_element_by_id('js-reactant-equivalent1').get_attribute("readonly")
        second_reactant_equivalent_field = self.driver.\
            find_element_by_id('js-reactant-equivalent2').get_attribute("readonly")
        self.assertFalse(first_reactant_equivalent_field)
        self.assertTrue(second_reactant_equivalent_field)
        """Now we check equivalents change with radiobutton"""
        first_reactant_equivalent = self.driver.\
            find_element_by_id('js-reactant-equivalent1').get_attribute("value")
        second_reactant_equivalent = self.driver.\
            find_element_by_id('js-reactant-equivalent2').get_attribute("value")
        first_reactant_equivalent_field = self.driver.\
            find_element_by_id('js-reactant-equivalent1').get_attribute("readonly")
        second_reactant_equivalent_field = self.driver.\
            find_element_by_id('js-reactant-equivalent2').get_attribute("readonly")
        self.assertEqual(first_reactant_equivalent, "0.5")
        self.assertEqual(second_reactant_equivalent, "1")
        self.assertFalse(first_reactant_equivalent_field)
        self.assertTrue(second_reactant_equivalent_field)
