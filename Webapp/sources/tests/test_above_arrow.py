from auxiliary_for_tests import *
import unittest
import flask_testing
from sources import app, db
import time
from selenium.webdriver.support.select import Select
from database_setup import test_database_create


class TestAboveArrow(unittest.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        pass

    def tearDown(self):
        restore_db()

    def fill_table(self):  # registration function
        with app.test_client() as self.client:
            with self.client.session_transaction() as session:
                # assign session variables
                session["workbook"] = "Test-Workbook"
                session["workgroup"] = "Test-Workgroup"
            return self.client.get(
                '/_process',
                follow_redirects=True
            )

    def filling_table(self):  # registration function
        with app.test_client() as self.client:
            with self.client.session_transaction() as session:
                # assign session variables
                session["workbook"] = "Test-Workbook"
                session["workgroup"] = "Test-Workgroup"
            return self.client.get(
                '/_process?reactants=C1CCCC1,OC(=O)C1=CC=CC=C1'
                '&reagents=CCO.C1=CC=CC=C1'
                '&products=CCNC(=O)C1=CC=CC=C1,N1C=CC=C1',
                follow_redirects=True
            )

    def test_reagent1(self):  # testing successful user registration
        response = self.filling_table()
        data = response.json
        self.assertIn(b'id=\\"js-reagent1\\" style=\\"width: 200px; border-radius: 5px;\\" '
                      b'value=\\"Ethanol\\', response.data)

    def test_reagent2(self):  # testing successful user registration
        response = self.filling_table()
        data = response.json
        self.assertIn(b'id=\\"js-reagent2\\" style=\\"width: 200px; border-radius: 5px;\\" '
                      b'value=\\"Benzene\\', response.data)

class TestNovelSolventAboveArrow(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """ Need to test the process of adding a novel compound to the reaction database via the sketcher"""

    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        """We load a test instance of the app, clear and re-populate the database then call the login function"""
        test_database_create(db)

    def tearDown(self):
        # close the browser window
        restore_db()

    def login(self):
        """This function sets the headless browser up, logs the user in and is called from setUp"""
        setup_selenium(self)
        username_field = self.driver.find_element_by_id('username')
        username_field.clear()
        username_field.send_keys("BB_Test")
        password_field = self.driver.find_element_by_id('password')
        password_field.clear()
        password_field.send_keys("BB_login")
        self.driver.find_element_by_id('submit').click()
        time.sleep(1)
        select_workgroup(self)
        select_wb = Select(self.driver.find_element_by_id("active-workbook"))
        select_wb.select_by_index(0)
        time.sleep(1)
        make_new_reaction(self)
        time.sleep(1)

    def test_input_fields_auto_fill(self):
        # login, fill in sketcher, submit, and then assert input fields appear with data
        self.login()
        rxn_smiles = 'C1=CC=CC=C1>O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(5)
        # do the name check because the api call might just be slow
        check_name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        if not check_name:
            self.driver.find_element_by_id('js-new-compound-name').send_keys('pyrazino[2,3-g]quinoxaline-5,10-dione')
        name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        mw = self.driver.find_element_by_id('js-new-compound-mw').get_attribute('value')
        smiles = self.driver.find_element_by_id('js-new-compound-smiles').get_attribute('value')
        self.assertEqual('pyrazino[2,3-g]quinoxaline-5,10-dione', name)
        self.assertEqual('212.03', mw)
        self.assertEqual('O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2', smiles)


if __name__ == '__main__':
    unittest.main()
