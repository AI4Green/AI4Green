from auxiliary_for_tests import *
import flask_testing
from sources import app, db
from database_setup import test_database_create


class HoverTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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

    # test hover
    def test_hover(self):
        hover_text = self.driver.find_element_by_class_name("js-add-reagent").get_attribute("title")
        self.assertIn(hover_text, "Add reagent to reaction")
        hover_text = self.driver.find_element_by_class_name("js-add-solvent").get_attribute("title")
        self.assertIn(hover_text, "Add solvent to reaction")
        hover_text = self.driver.find_element_by_id("js-add-new-reagent-by-table").get_attribute("title")
        self.assertIn(hover_text, "Add new reagent to database")
