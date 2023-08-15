from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# multiprocessing.set_start_method("fork")

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

    # test reaction table displays solvent
    def test_solvent_not_there_initially(self):
        pg_source = self.driver.page_source
        self.assertNotIn('js-solvent1', pg_source)

    # test reaction table displays solvent but not solvent until added
    def test_solvent_in_reaction_table(self):
        demo_reaction(self)
        time.sleep(3)
        pg_source = self.driver.page_source
        self.assertNotIn('js-solvent1', pg_source)
        self.assertIn('Solvents', pg_source)

    # test solvent buttons works
    def test_solvent_buttons(self):
        demo_reaction(self)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('js-solvent1', pg_source)
        self.assertNotIn('js-solvent2', pg_source)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('js-solvent1', pg_source)
        self.assertIn('js-solvent2', pg_source)
        self.driver.find_element_by_class_name('js-remove-solvent').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('js-solvent1', pg_source)
        self.assertNotIn('js-solvent2', pg_source)

    # test select solvent works
    def test_select_solvent(self):
        demo_reaction(self)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # test background colour is not initially grey
        rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
        self.assertNotEqual(rgb, "rgba(139, 0, 0, 1)")
        self.driver.find_element_by_id('js-solvent1').click()
        self.driver.find_element_by_id('js-solvent1').send_keys("ben")
        time.sleep(1)
        move_down_one_option_action_chains(self)
        time.sleep(1)
        # after "ben" typed the first entry should have value "Benzene" and have grey background
        value = self.driver.find_element_by_id('js-solvent1').get_attribute("value")
        self.assertEqual(value, "Benzene")
        rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
        self.assertEqual(rgb, "rgba(139, 0, 0, 1)")

    # test by selecting second solvent
    def test_select_solvent2(self):
        demo_reaction(self)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # test background colour is not initially green
        rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
        self.assertNotEqual(rgb, "rgba(0, 255, 0, 1))")
        self.driver.find_element_by_id('js-solvent1').click()
        self.driver.find_element_by_id('js-solvent1').send_keys("ben")
        time.sleep(1)
        actions = webdriver.ActionChains(self.driver)
        actions.send_keys(Keys.ARROW_DOWN)
        actions.send_keys(Keys.ARROW_DOWN)
        actions.send_keys(Keys.ENTER)
        actions.perform()
        time.sleep(1)
        # after "ben" typed the second entry should have value "Benzyl alcohol" and have green background
        value = self.driver.find_element_by_id('js-solvent1').get_attribute("value")
        self.assertEqual(value, "Benzyl alcohol")
        rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
        self.assertEqual(rgb, "rgba(0, 255, 0, 1)")

    def test_second_solvent(self):
        """
        This tests second solvent can be added and the first solvent can be updated without affecting the other solvent
        """
        demo_reaction(self)
        time.sleep(3)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # Adding the first solvent
        rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
        self.assertNotEqual(rgb, "rgba(139, 0, 0, 1)")
        self.driver.find_element_by_id('js-solvent1').click()
        self.driver.find_element_by_id('js-solvent1').send_keys("Benzene")
        time.sleep(1)
        move_down_one_option_action_chains(self)
        time.sleep(1)
        # after "ben" typed the first entry should have value "Benzene" and have grey background
        value = self.driver.find_element_by_id('js-solvent1').get_attribute("value")
        self.assertEqual(value, "Benzene")
        rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
        time.sleep(1)
        self.assertEqual(rgb, "rgba(139, 0, 0, 1)")
        # Adding the second solvent
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(2)
        rgb = self.driver.find_element_by_id('js-solvent2').value_of_css_property('background-color')
        self.assertNotEqual(rgb, "rgba(0, 255, 0, 1))")
        self.driver.find_element_by_id('js-solvent2').click()
        self.driver.find_element_by_id('js-solvent2').send_keys("Benzyl alcohol")
        move_down_one_option_action_chains(self)
        time.sleep(2)
        # the second entry should have value "Benzyl alcohol" and have green background
        value = self.driver.find_element_by_id('js-solvent2').get_attribute("value")
        self.assertEqual(value, "Benzyl alcohol")
        rgb = self.driver.find_element_by_id('js-solvent2').value_of_css_property('background-color')
        self.assertEqual(rgb, "rgba(0, 255, 0, 1)")
        # Now change the first solvent to dichloromethane, a red coded solvent
        self.driver.find_element_by_id('js-solvent1').click()
        self.driver.find_element_by_id('js-solvent1').clear()
        self.driver.find_element_by_id('js-solvent1').send_keys("Dichloromethane")
        move_down_one_option_action_chains(self)
        # now check js-solvent2 hasn't been changed by updating js-solvent1
        value = self.driver.find_element_by_id('js-solvent2').get_attribute("value")
        self.assertEqual(value, "Benzyl alcohol")
        rgb = self.driver.find_element_by_id('js-solvent2').value_of_css_property('background-color')
        self.assertEqual(rgb, "rgba(0, 255, 0, 1)")
        # check the newly added solvent field has the expected properties
        time.sleep(2)
        value = self.driver.find_element_by_id('js-solvent1').get_attribute("value")
        self.assertEqual(value, "Dichloromethane")
        rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
        self.assertEqual(rgb, "rgba(255, 0, 0, 1)")

    # test new solvent cas works
    def test_new_solvent_cas(self):
        demo_reaction(self)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        self.driver.find_element_by_id('js-solvent1').click()
        self.driver.find_element_by_id('js-solvent1').send_keys("123-45-6")
        time.sleep(1)
        self.driver.switch_to.alert.accept()
        time.sleep(1)
        self.driver.find_element_by_id('js-solvent-rounded-concentration1').click()
        time.sleep(1)
        # after CAS typed it should in the CAS field for a new reagent
        value = self.driver.find_element_by_id('js-input-solvent-cas').get_attribute("value")
        self.assertEqual(value, "123-45-6")

