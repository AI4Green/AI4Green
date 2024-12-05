from auxiliary_for_tests import *
from selenium.webdriver.common.keys import Keys
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create


# test solvent selection
class SolventSelectionTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    # test reaction table displays reagent
    def test_reagent_not_there_initially(self):
        pg_source = self.driver.page_source
        self.assertNotIn('Catalysts/reagents', pg_source)

    # test reaction table displays solvent but not solvent until added
    def test_reagent_in_reaction_table(self):
        demo_reaction(self)
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertNotIn('js-reagent1', pg_source)
        self.assertIn('Catalysts/reagents', pg_source)

    # test reagent appears only in its input field
    def test_reagent_input(self):
        demo_reaction(self)
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_id('js-reagent1').send_keys('Aniline')
        self.driver.find_element_by_id('js-reagent1').send_keys(Keys.ENTER)
        reagent1 = self.driver.find_element_by_id('js-reagent1').get_attribute('value')
        time.sleep(1)
        reagent2 = self.driver.find_element_by_id('js-reagent2').get_attribute('value')
        self.assertNotEqual('Aniline', reagent2)

    # test add/remove reagent buttons works
    def test_reagent_buttons(self):
        demo_reaction(self)
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('js-reagent1', pg_source)
        self.assertNotIn('js-reagent2', pg_source)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('js-reagent1', pg_source)
        self.assertIn('js-reagent2', pg_source)
        self.driver.find_element_by_class_name('js-remove-reagent').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('js-reagent1', pg_source)
        self.assertNotIn('js-reagent2', pg_source)

    # test reagent data returned
    def test_reagent_data_returned(self):
        demo_reaction(self)
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        reagent = self.driver.find_element_by_id('js-reagent1')
        reagent.clear()
        reagent.send_keys('Aniline')
        reagent_molecular_weight = self.driver.find_element_by_id('js-reagent-molecular-weight1')
        reagent_molecular_weight.click()
        time.sleep(1)
        reagent_density = self.driver.find_element_by_id('js-reagent-density1')
        reagent_concentration = self.driver.find_element_by_id('js-reagent-concentration1')
        reagent_hazards = self.driver.find_element_by_id('js-reagent-hazards1')
        self.assertIn("93.13", reagent_molecular_weight.get_attribute('value'))
        self.assertIn("1.022", reagent_density.get_attribute('value'))
        self.assertIn("H301-H311-H317-H318-H331-H341-H351-H372-H400",
                      reagent_hazards.get_attribute('value'))
