from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time
import os
from utilities import basedir
import unittest


class TestMarvinPresent(flask_testing.TestCase):
    """These functions test the Marvin sketch tool and key for use are present in the response from the webpage."""

    def login(self, username, password):  # login function
        return self.client.post(
            '/auth/login',
            data=dict(username=username, password=password),
            follow_redirects=True
        )

    def logout(self):  # logout function
        return self.client.get(
            '/auth/logout',
            follow_redirects=True
        )

    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        """We login each time"""
        self.response = self.login('PI_Test', 'PI_Login')

    def tearDown(self):
        self.logout()
        restore_db()

    def test_marvin_key(self):
        """
        Here we test we have the marvin key in the webpage response and this matches the key in the base html file.
        """
        base_html_file = os.path.join(basedir, 'sources/templates/base.html')
        with open(base_html_file, 'r') as f:
            html_string = f.read()
        # old api key
        # self.assertIn(b'https://marvinjs.chemicalize.com/v1/69b0c6c7851547258e92db7880f174ef/client-settings.js',
        #               self.response.data)
        # self.assertIn('https://marvinjs.chemicalize.com/v1/69b0c6c7851547258e92db7880f174ef/client-settings.js',
        #               html_string)
        self.assertIn(b'"https://marvinjs.chemicalize.com/v1/bb45b77a07a64428b9825de4be63c424/client-settings.js"',
                      self.response.data)
        self.assertIn('"https://marvinjs.chemicalize.com/v1/{{ marvin_js_key }}/client-settings.js"',
                      html_string)

    def test_marvin_window(self):
        self.assertIn(b'https://marvinjs.chemicalize.com/v1/client.js', self.response.data)


class MarvinInteractionTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """
    The next tests use selenium and the LiveServerTestCase to open a browser and check Marvin js behaves as expected
    """

    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        """We load a test instance of the app, clear the db register a new user then login"""
        setup_selenium(self)
        login(self)
        self.driver.find_element_by_id("TopNavSketcherButton").click()
        time.sleep(2)

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def test_marvin_loads(self):
        """This function tests that the marvin window is present in the page source"""
        pg_source = self.driver.page_source
        self.assertIn('marvin-test', pg_source, "The marvin window must be present")

    def test_marvin_demo_function(self):
        """Here we test that compounds can be loaded into marvin_js by the example demo button by pressing the submit
         button and checking the correct compound appears in the reaction table"""
        pg_source = self.driver.page_source
        self.assertNotIn('Benzoic acid', pg_source, "no compound names should be present at this stage")
        time.sleep(3)
        demo_reaction(self)
        pg_source = self.driver.page_source
        self.assertIn('Benzoic acid', pg_source, "benzoic acid should now be loaded in the reaction table")
        reactant_mass_field = self.driver.find_element_by_id("js-reactant-rounded-mass1")
        reactant_mass_field.send_keys(100)
        self.driver.find_element_by_id('js-reactant-equivalent2').send_keys(1)
        phys_form_dic = {'js-reactant-physical-form1': 1, 'js-reactant-physical-form2': 1,
                         'js-product-physical-form1': 1}
        fill_in_physical_forms(self, phys_form_dic)
        time.sleep(2)
        self.driver.find_element_by_id('action-summary').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Sustainability (CHEM21)', pg_source)


if __name__ == '__main__':
    unittest.main()
