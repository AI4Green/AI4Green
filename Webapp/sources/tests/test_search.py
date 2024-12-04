from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create


# test reaction search
class SearchTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self, headless='no')
        login(self, "BB_Test", "BB_login")

    def test_exact_structure_search(self):
        self.driver.find_element_by_id('TopNavSearchButton').click()
        time.sleep(1)
        add_reaction_sketcher(self, 'O=N(=O)C1=CC=CC=C1')
        self.driver.find_element_by_id('structure-search-btn').click()
        time.sleep(2)
        search_results = self.driver.find_element_by_id('search-results-message').get_attribute('innerHTML')
        self.assertIn('2 results found', search_results)
        reaction_name = self.driver.find_element_by_id('reaction-name1').text
        reaction_name2 = self.driver.find_element_by_id('reaction-name2').text
        self.assertEqual(['a nitro reduction2', 'nitro reduction2'], [reaction_name, reaction_name2])

    def test_reactions_reload(self):
        self.driver.find_element_by_id('TopNavSearchButton').click()
        time.sleep(1)
        add_reaction_sketcher(self, 'O=N(=O)C1=CC=CC=C1')
        self.driver.find_element_by_id('structure-search-btn').click()
        time.sleep(2)
        self.driver.find_element_by_id('a nitro reduction2-reload').click()
        time.sleep(5)
        # confirm current tab is still search page
        search_results = self.driver.find_element_by_id('search-results-message').get_attribute('innerHTML')
        self.assertIn('2 results found', search_results)
        # now change tab and confirm reaction page
        new_window = self.driver.window_handles[1]
        self.driver.switch_to_window(new_window)
        reaction_name = self.driver.find_element_by_id('js-reaction-name2').get_attribute('value')
        self.assertEqual('a nitro reduction2', reaction_name)
