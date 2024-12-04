from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time


# test solvent guide
class SolventGuideTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    def test_top_banner_button(self):
        """Tests the solvent guide page renders when the top banner button is pressed"""
        self.driver.find_element_by_id("solvent-guide").click()
        self.assertIn("<h2>Solvent Guide</h2>", self.driver.page_source)

    def test_correct_cards_renders(self):
        """Tests the correct card is active when a solvent is selected"""
        self.driver.find_element_by_id("solvent-guide").click()
        time.sleep(2)
        self.driver.find_element_by_id("family1").click()
        self.driver.find_element_by_id("solvent5").click()
        time.sleep(1)
        self.assertIn('active show" id="list-50"', self.driver.page_source)
        self.driver.find_element_by_id("family2").click()
        self.driver.find_element_by_id("solvent11").click()
        time.sleep(1)
        self.assertIn('active show" id="list2-7"', self.driver.page_source)
        self.driver.find_element_by_id("family5").click()
        self.driver.find_element_by_id("solvent24").click()
        time.sleep(1)
        self.assertIn('active show" id="list-14"', self.driver.page_source)
        self.assertNotIn('active show" id="list-50"', self.driver.page_source)

    def test_suggestion_opens_opposite_card(self):
        """Tests suggested solvent opens on opposite card"""
        self.driver.find_element_by_id("solvent-guide").click()
        time.sleep(2)
        self.driver.find_element_by_id("family7").click()
        time.sleep(1)
        self.driver.find_element_by_id("solvent15").click()
        time.sleep(1)
        self.driver.find_element_by_id("sub1_15").click()
        time.sleep(1)
        self.assertIn('active show" id="list2-33"', self.driver.page_source)
        self.assertIn('active show" id="list-34"', self.driver.page_source)

    def test_preload_no_solvent_reaction_table(self):
        """Tests no solvent is preloaded if solvent in CHEM21 not selected before pressing solvent guide button in reaction table"""
        select_workgroup(self)
        select_workbook(self)
        make_new_reaction(self)
        demo_reaction(self)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        self.driver.find_element_by_id('go-to-solvent-guide1').click()
        time.sleep(2)
        # switch to new tab
        self.driver.switch_to.window(self.driver.window_handles[1])
        self.assertIn("<h2>Solvent Guide</h2>", self.driver.page_source)
        # test no solvents loaded initally
        self.assertIn('active show" id="list-default"', self.driver.page_source)

    def test_preload_solvent_reaction_table(self):
        """Tests solvent is preloaded if solvent in CHEM21 selected before pressing solvent guide button in reaction table"""
        select_workgroup(self)
        select_workbook(self)
        make_new_reaction(self)
        demo_reaction(self)
        self.driver.find_element_by_class_name('js-add-solvent').click()
        time.sleep(1)
        # get benzene in solvent box
        self.driver.find_element_by_id('js-solvent1').click()
        self.driver.find_element_by_id('js-solvent1').send_keys("ben")
        time.sleep(1)
        move_down_one_option_action_chains(self)

        # actions = webdriver.ActionChains(self.driver)
        # actions.send_keys(Keys.ARROW_DOWN)
        # actions.send_keys(Keys.ENTER)
        # actions.perform()
        time.sleep(1)
        self.driver.find_element_by_id('go-to-solvent-guide1').click()
        time.sleep(2)
        # switch to new tab
        self.driver.switch_to.window(self.driver.window_handles[1])
        self.assertIn("<h2>Solvent Guide</h2>", self.driver.page_source)
        # test benzene loaded initially
        self.assertIn('active show" id="list-30"', self.driver.page_source)
