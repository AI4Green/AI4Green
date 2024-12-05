from auxiliary_for_tests import *
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import unittest
import flask_testing
import multiprocessing
from database_setup import test_database_create
from sources import app, db
import time



#multiprocessing.set_start_method("fork")

def register(self, fullname, username, email, password, password2):  # registration function
    return self.app.post(
        '/auth/register',
        data=dict(fullname=fullname, username=username, email=email, password=password, password2=password2),
        follow_redirects=True
    )


class LoadSummaryTableTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):

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
        demo_reaction(self)
        # we need to fill in the mass of the limiting reactant
        reactant_mass_field = self.driver.find_element_by_id("js-reactant-rounded-mass1")
        reactant_mass_field.send_keys(100)
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys(1)

    def test_summary_function(self):
        """Here we test that compounds can be loaded into marvin_js by the example demo button by pressing the submit
         button and checking the correct compound appears in the reaction table"""

        # action summary button loads the summary table
        phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                         'js-product-physical-form1': 1}
        fill_in_physical_forms(self, phys_form_dic)
        self.driver.find_element_by_id("action-summary").click()
        time.sleep(3)
        pg_source = self.driver.page_source

        self.assertIn('Theoretical Yield' and 'Sustainability', pg_source)

        # Fill in yield
        yield_field = self.driver.find_element_by_id('js-real-product-mass')
        yield_field.send_keys(400)
        yield_str = 'js-yield" style="width: 80px; border: none; background-color: rgb(0, 255, 0);'
        yield_field.send_keys(Keys.TAB)
        pg_source = self.driver.page_source
        self.assertIn(yield_str, pg_source)

        # Fill in temperature
        temp_field = self.driver.find_element_by_id('js-temperature')
        temp_field.send_keys(100)
        temp_str = 'id="js-temperature-cell" class="hazard-warning" style="background-color: rgb(255, 255, 0);'
        temp_field.send_keys(Keys.TAB)
        pg_source = self.driver.page_source
        self.assertIn(temp_str, pg_source)

        # sustainability dropdown for isolation
        dropdown1 = self.driver.find_element_by_id('js-isolation')
        dropdown1.click()
        actions = webdriver.ActionChains(self.driver)
        actions.send_keys(Keys.ARROW_DOWN)
        actions.send_keys(Keys.ENTER)
        actions.perform()
        time.sleep(1)
        pg_source = self.driver.page_source

        # Test Hazard risk-score function
        self.driver.find_element_by_id('major').click()
        self.driver.find_element_by_id('frequentOccur').click()
        self.driver.find_element_by_id('buildingWide').click()
        self.assertTrue(self.driver.find_element_by_id('categoryA').is_selected())

        # Assert is below to give time for style change to red to appear in HTML pg source
        self.assertIn('id="js-isolation-cell" class="hazard-hazardous" style="background-color: rgb(255, 0, 0);',
                      pg_source)

        # Test other hazards text box can be written in
        try:
            risks_textbox = self.driver.find_element_by_id('other-risks-textbox')
            risks_textbox.send_keys('Lachrymator used, rinse glassware thoroughly before '
                                    'removing from fumehood')
        except:
            self.fail("writing in the other risks textbox raised an exception unexpectedly!")
        # Confirms print button element is present
        self.assertTrue(self.driver.find_element_by_id('print-pdf'))


if __name__ == '__main__':
    unittest.main()
