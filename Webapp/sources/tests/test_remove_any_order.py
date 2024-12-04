from auxiliary_for_tests import *
import flask_testing
from sources import app, db
from database_setup import test_database_create


# test remove any order
class RemoveAnyOrderTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        demo_reaction(self)

    # test remove reagent any order
    def test_remove_any_order(self):
        # add three reagents and three solvents
        self.driver.find_elements_by_class_name('js-add-reagent')[0].click()
        self.driver.find_elements_by_class_name('js-add-reagent')[0].click()
        self.driver.find_elements_by_class_name('js-add-reagent')[0].click()
        self.driver.find_elements_by_class_name('js-add-solvent')[0].click()
        self.driver.find_elements_by_class_name('js-add-solvent')[0].click()
        self.driver.find_elements_by_class_name('js-add-solvent')[0].click()
        # add some equivalence to these reagents
        self.driver.find_element_by_id('js-reagent-equivalent1').send_keys('100')
        self.driver.find_element_by_id('js-reagent-equivalent2').send_keys('200')
        self.driver.find_element_by_id('js-reagent-equivalent3').send_keys('300')
        # add some volumes to these solvents
        self.driver.find_element_by_id('js-solvent-volume1').send_keys('10')
        self.driver.find_element_by_id('js-solvent-volume2').send_keys('20')
        self.driver.find_element_by_id('js-solvent-volume3').send_keys('30')
        # delete middle reagent
        self.driver.find_element_by_id('remove-reagent2').click()
        # make sure first and third reagent are still there
        eq1 = self.driver.find_element_by_id('js-reagent-equivalent1').get_attribute('value')
        eq2 = self.driver.find_element_by_id('js-reagent-equivalent2').get_attribute('value')
        self.assertEqual(eq1, "100")
        self.assertEqual(eq2, "300")
        # check id and table numbers are updated
        table_number = self.driver.find_element_by_id('js-reagent-table-number2').get_attribute('value')
        self.assertEqual(table_number, "4")
        table_number = self.driver.find_element_by_id('js-solvent-table-number1').get_attribute('value')
        self.assertEqual(table_number, "5")
        table_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute('value')
        self.assertEqual(table_number, "8")
        # remove first solvent
        self.driver.find_element_by_id('remove-solvent1').click()
        # make sure first and third reagent are still there
        eq1 = self.driver.find_element_by_id('js-solvent-volume1').get_attribute('value')
        eq2 = self.driver.find_element_by_id('js-solvent-volume2').get_attribute('value')
        self.assertEqual(eq1, "20")
        self.assertEqual(eq2, "30")
        # check id and table numbers are updated
        table_number = self.driver.find_element_by_id('js-solvent-table-number2').get_attribute('value')
        self.assertEqual(table_number, "6")
        table_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute('value')
        self.assertEqual(table_number, "7")
