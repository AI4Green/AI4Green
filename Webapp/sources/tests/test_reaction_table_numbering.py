from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time

# Test Add/Remove from WB/WG on the manage workbook page
class ReactionNumberingTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """This test uses selenium to test add/remove feature on manage workgroup page"""

    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        setup_selenium(self)
        login(self)
        select_workgroup(self)
        select_workbook(self)
        make_new_reaction(self)
        demo_reaction(self)

    def tearDown(self):
        restore_db()
        self.driver.quit()

    def test_initial_numbering(self):
        """Tests the initial numbering"""
        # product initially has value = 3 and 4
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        # p2_number = self.driver.find_element_by_id('js-product-table-number2').get_attribute("value")
        self.assertEqual(p1_number, "3")
        # self.assertEqual(p2_number, "4")
        # nothing has value 5
        self.assertNotIn("""value=\"4\"""", self.driver.page_source)

    def test_add_reagent(self):
        """Tests adding a reagent"""
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        # new reagent has value = 3
        r1_number = self.driver.find_element_by_id('js-reagent-table-number1').get_attribute("value")
        self.assertEqual(r1_number, "3")
        # product now has value = 4
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "4")
        # nothing has value 6
        self.assertNotIn("""value=\"6\"""", self.driver.page_source)
        """Tests adding a second reagent"""
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        # first reagent has value = 3
        r1_number = self.driver.find_element_by_id('js-reagent-table-number1').get_attribute("value")
        self.assertEqual(r1_number, "3")
        # second reagent has value = 4
        r2_number = self.driver.find_element_by_id('js-reagent-table-number2').get_attribute("value")
        self.assertEqual(r2_number, "4")
        # product now has value = 5
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "5")
        # nothing has value 7
        self.assertNotIn("""value=\"7\"""", self.driver.page_source)

    def test_remove_reagent(self):
        """Tests removing reagents"""
        # add two reagents
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        # remove reagent and check values
        #self.driver.find_element_by_class_name('js-remove-reagent').click()
        self.driver.find_element_by_id('remove-reagent1').click()
        time.sleep(1)
        # reagent has value = 3
        r1_number = self.driver.find_element_by_id('js-reagent-table-number1').get_attribute("value")
        self.assertEqual(r1_number, "3")
        # product now has value = 4
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "4")
        # nothing has value 6
        self.assertNotIn("""value=\"6\"""", self.driver.page_source)
        """Tests removing a second reagent"""
        self.driver.find_element_by_id('remove-reagent1').click()
        time.sleep(1)
        # product now has value = 3
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "3")
        # nothing has value 5
        self.assertNotIn("""value=\"5\"""", self.driver.page_source)

    def test_add_solvent(self):
        """Tests adding a solvent"""
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # new solvent has value = 2
        s1_number = self.driver.find_element_by_id('js-solvent-table-number1').get_attribute("value")
        self.assertEqual(s1_number, "3")
        # product now has value = 3
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "4")
        # nothing has value 4
        self.assertNotIn("""value=\"6\"""", self.driver.page_source)
        """Tests adding a second solvent"""
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # first solvent has value = 3
        s1_number = self.driver.find_element_by_id('js-solvent-table-number1').get_attribute("value")
        self.assertEqual(s1_number, "3")
        # second solvent has value = 4
        s2_number = self.driver.find_element_by_id('js-solvent-table-number2').get_attribute("value")
        self.assertEqual(s2_number, "4")
        # product now has value = 5
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "5")
        # nothing has value 7
        self.assertNotIn("""value=\"7\"""", self.driver.page_source)

    def test_remove_solvent(self):
        """Tests removing solvents"""
        # add two reagents
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # remove solvent and check values
        self.driver.find_element_by_class_name('js-remove-solvent').click()
        time.sleep(1)
        # solvent has value = 3
        s1_number = self.driver.find_element_by_id('js-solvent-table-number1').get_attribute("value")
        self.assertEqual(s1_number, "3")
        # product now has value = 4
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "4")
        # nothing has value 6
        self.assertNotIn("""value=\"6\"""", self.driver.page_source)
        """Tests removing a second solvent"""
        self.driver.find_element_by_class_name('js-remove-solvent').click()
        time.sleep(1)
        # product now has value = 3
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "3")
        # nothing has value 5
        self.assertNotIn("""value=\"5\"""", self.driver.page_source)

    def test_mixture_add_remove_reagents_solvents(self):
        """Tests the response of adding or removing reagents or solvents"""
        # add two reagents and two solvents
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # check values
        r1_number = self.driver.find_element_by_id('js-reagent-table-number1').get_attribute("value")
        self.assertEqual(r1_number, "3")
        r2_number = self.driver.find_element_by_id('js-reagent-table-number2').get_attribute("value")
        self.assertEqual(r2_number, "4")
        s1_number = self.driver.find_element_by_id('js-solvent-table-number1').get_attribute("value")
        self.assertEqual(s1_number, "5")
        s2_number = self.driver.find_element_by_id('js-solvent-table-number2').get_attribute("value")
        self.assertEqual(s2_number, "6")
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "7")
        # remove one reagent and add one solvent
        self.driver.find_element_by_class_name('js-remove-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # check values
        r1_number = self.driver.find_element_by_id('js-reagent-table-number1').get_attribute("value")
        self.assertEqual(r1_number, "3")
        s1_number = self.driver.find_element_by_id('js-solvent-table-number1').get_attribute("value")
        self.assertEqual(s1_number, "4")
        s2_number = self.driver.find_element_by_id('js-solvent-table-number2').get_attribute("value")
        self.assertEqual(s2_number, "5")
        s3_number = self.driver.find_element_by_id('js-solvent-table-number3').get_attribute("value")
        self.assertEqual(s3_number, "6")
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "7")
        # add two reagents and remove one solvent
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-add-reagent').click()
        time.sleep(1)
        self.driver.find_element_by_class_name('js-remove-solvent').click()
        time.sleep(1)
        # check values
        r1_number = self.driver.find_element_by_id('js-reagent-table-number1').get_attribute("value")
        self.assertEqual(r1_number, "3")
        r2_number = self.driver.find_element_by_id('js-reagent-table-number2').get_attribute("value")
        self.assertEqual(r2_number, "4")
        r3_number = self.driver.find_element_by_id('js-reagent-table-number3').get_attribute("value")
        self.assertEqual(r3_number, "5")
        s1_number = self.driver.find_element_by_id('js-solvent-table-number1').get_attribute("value")
        self.assertEqual(s1_number, "6")
        s2_number = self.driver.find_element_by_id('js-solvent-table-number2').get_attribute("value")
        self.assertEqual(s2_number, "7")
        p1_number = self.driver.find_element_by_id('js-product-table-number1').get_attribute("value")
        self.assertEqual(p1_number, "8")
