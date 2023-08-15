from auxiliary_for_tests import *
from pony.orm import db_session, select
import flask_testing
from database_setup import test_database_create
from sources import app, db
import time
from selenium.webdriver.support.select import Select


class CompoundDataErrorReportTestCase(flask_testing.TestCase):
    """
    Tests compound error reports for incorrect data.
    """

    def create_app(self):
        app.config.from_object('config.TestConfig')
        app.config['LOGIN_DISABLED'] = True
        test_database_create(db)
        return app

    def setUp(self):
        """We load a test instance of the app, clear and re-populate the database then call the login function"""
        pass

    def tearDown(self):
        restore_db()

    def report_error(self, compound_name, error_type, additional_info, compound_ID):
        """Posts data to the compound_data_error routes"""
        return self.client.post(
            '/compound_data_error_report',
            data=dict(compoundName=compound_name, errorType=error_type, additionalInfo=additional_info, compoundID=compound_ID)
        )

    @db_session
    def test_error_reported(self):
        """Test an error is stored in the db"""
        compound = select(x for x in db.Compound if x.id == 68).first()
        response = self.report_error("Benzoic acid", "incorrect-cas", "unit test", "68")
        error_details = select((x.compound_name, x.error_type, x.additional_info, x.compound) for x in db.CompoundDataErrorReport if x.additional_info == "unit test").first()
        self.assertEqual(("Benzoic acid", "incorrect-cas", "unit test", compound), error_details)


class ReportCompoundDataErrorLiveTestCase(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """This tests the report compound data functionality"""
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        setup_selenium(self)
        login(self, "BB_Test", "BB_login")
        select_workgroup(self)

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def make_reaction_and_report_reactant1(self):
        make_new_reaction(self)
        demo_reaction(self)
        self.driver.find_element_by_id("js-report-reactant1").click()
        time.sleep(1)

    def test_compound_name(self):
        self.make_reaction_and_report_reactant1()
        title = self.driver.find_element_by_id("compound-data-report-modal-title").text
        self.assertEqual('Benzoic acid', title)

    def test_dropdown(self):
        self.make_reaction_and_report_reactant1()
        dropdown = Select(self.driver.find_element_by_id("compound-data-report-dropdown")).select_by_visible_text("Incorrect molecular weight")
        dropdown_option = self.driver.find_element_by_id("compound-data-report-dropdown").get_attribute('value')
        self.assertEqual('incorrect-molecular-weight', dropdown_option)

    def test_textbox(self):
        self.make_reaction_and_report_reactant1()
        self.driver.find_element_by_id('compound-data-report-additional-info').send_keys("Weight is too low")
        text_box = self.driver.find_element_by_id('compound-data-report-additional-info').get_attribute('value')
        self.assertEqual('Weight is too low', text_box)

    @db_session
    def test_save_report(self):
        self.make_reaction_and_report_reactant1()
        self.driver.find_element_by_id('compound-data-report-additional-info').send_keys("unit test")
        dropdown = Select(self.driver.find_element_by_id("compound-data-report-dropdown")).select_by_visible_text(
            "Incorrect CAS")
        self.driver.find_element_by_id("report-compound-data-submit").click()
        time.sleep(1)
        submit_message = self.driver.find_element_by_id('compound-data-report-text').text
        self.assertEqual('Thank you for your report', submit_message)
        compound = select(x for x in db.Compound if x.name.lower() == 'benzoic acid').first()
        error_details = select((x.compound_name, x.error_type, x.additional_info, x.compound) for x in db.CompoundDataErrorReport if x.additional_info == "unit test").first()
        self.assertEqual(("Benzoic Acid", "incorrect-cas", "unit test", compound), error_details)

    def test_buttons_other_components(self):
        """This tests that the buttons appear and work for solvents/reagents/products and they clear between uses"""
        self.driver.find_element_by_id("nitro reduction2-reload").click()
        time.sleep(3)
        self.driver.find_element_by_id("js-report-reactant2").click()
        time.sleep(1)
        self.driver.find_element_by_id('compound-data-report-additional-info').send_keys("unit test")
        dropdown = Select(self.driver.find_element_by_id("compound-data-report-dropdown")).select_by_visible_text(
            "Incorrect molecular weight")
        self.driver.find_element_by_id("report-compound-data-close").click()
        time.sleep(1)
        self.driver.find_element_by_id("js-report-product1").click()
        time.sleep(1)
        text_box = self.driver.find_element_by_id('compound-data-report-additional-info').get_attribute('value')
        dropdown_option = self.driver.find_element_by_id("compound-data-report-dropdown").get_attribute('value')
        self.assertEqual('', text_box)
        self.assertEqual('-select-', dropdown_option)
        self.driver.find_element_by_id("report-compound-data-submit").click()
        text_box = self.driver.find_element_by_id('compound-data-report-additional-info').get_attribute('value')
        dropdown_option = self.driver.find_element_by_id("compound-data-report-dropdown").get_attribute('value')
        self.assertEqual('', text_box)
        self.assertEqual('-select-', dropdown_option)

