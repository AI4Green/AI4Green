from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create


# data loss warnings when resubmitting reaction data
class DataLossWarningTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    def test_resubmit_warning(self):
        """Tests that resubmitting the reaction or summary table gives a data loss warning"""
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
        # resubmit, dismiss warning, and check data is still there
        self.driver.execute_script("arguments[0].scrollIntoView();", self.driver.find_element_by_id('demo-button'))
        time.sleep(1)
        self.driver.find_element_by_id('action-button-submit').click()
        self.driver.switch_to.alert.dismiss()
        self.assertEqual(self.driver.find_element_by_id('js-reactant-rounded-mass1').get_attribute("value"), "500")
        # resubmit, check warning and accept warning, check data gone
        self.driver.find_element_by_id('action-button-submit').click()
        self.driver.switch_to.alert.accept()
        time.sleep(2)
        self.assertEqual(self.driver.find_element_by_id('js-reactant-rounded-mass1').get_attribute("value"), "")
        # re-add in information
        self.driver.find_element_by_id('js-reactant-rounded-mass1').send_keys("500")
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys("2")
        phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                         'js-product-physical-form1': 1}
        fill_in_physical_forms(self, phys_form_dic)
        # press summary button
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(2)
        # add some information
        self.driver.find_element_by_id('js-temperature').send_keys("80")
        # resubmit, dismiss warning, and check data is still there
        self.driver.execute_script("arguments[0].scrollIntoView();", self.driver.find_element_by_id('js-add-new-reagent-by-table'))
        time.sleep(1)
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        self.driver.switch_to.alert.dismiss()
        self.assertEqual(self.driver.find_element_by_id('js-temperature').get_attribute("value"), "80")
        # resubmit, check warning and accept warning, check data gone
        self.driver.find_element_by_id('action-summary').click()
        self.driver.switch_to.alert.accept()
        time.sleep(2)
        self.assertEqual(self.driver.find_element_by_id('js-temperature').get_attribute("value"), "")

