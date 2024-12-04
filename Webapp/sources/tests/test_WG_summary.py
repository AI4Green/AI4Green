from auxiliary_for_tests import *
from selenium.webdriver.support.select import Select
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create
import unittest


# function to log in
def login_routes(self, username, password):
    return self.client.post(
        '/auth/login',
        data=dict(username=username, password=password),
        follow_redirects=True
    )


# function to log out
def logout_routes(self):
    return self.client.get(
        '/auth/logout',
        follow_redirects=True
    )

def WG_summary_page(self):
    return self.client.post(
        '/workgroup_membership_summary',
    follow_redirects=True
    )


# test WG summary page can only be accessed when logged in
class TestReturnToLogin(flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        app.config['LOGIN_DISABLED'] = False
        test_database_create(db)
        return app

    def setUp(self):
        login_routes(self, 'SR_Test', 'SR_login')

    def tearDown(self):
        restore_db()

    # test we return to log in page if not logged in
    def test_not_logged_in(self):
        logout_routes(self)
        response_WG_summary_page = WG_summary_page(self)
        self.assertIn(b'Please log in to access this page', response_WG_summary_page.data)

    # test we get to the WG summary page if logged in
    def test_get_to_WG_summary(self):
        response_WG_summary_page = WG_summary_page(self)
        self.assertIn(b'Workgroup Membership Summary', response_WG_summary_page.data)


# test correct information is shown when selecting and displaying workgroups/workbooks
class WGInteractionTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        self.driver.find_element_by_id('WG-membership-summary').click()
        time.sleep(1)

    # test summary page loads
    def test_WG_summary_loads(self):
        pg_source = self.driver.page_source
        self.assertIn('Workgroup Membership Summary', pg_source)

    # test that new text does not appear before you select a workgroup and press button
    def test_before_select(self):
        self.driver.find_element_by_id('submit').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertNotIn('You are a member of the following workbooks', pg_source)

    # test that new text appears when you select a workgroup and press button
    def test_select_option(self):
        select_wg = Select(self.driver.find_element_by_id("workgroups"))
        select_wg.select_by_index(1)
        self.driver.find_element_by_id('submit').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('You are a member of the following workbooks', pg_source)

    # test that correct information is shown when a workgroup is selected
    def test_correct_information(self):
        select_wg = Select(self.driver.find_element_by_id("workgroups"))
        select_wg.select_by_index(1)
        self.driver.find_element_by_id('submit').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertNotIn('Test-Workbook3', pg_source)
        self.assertIn('Test-Workbook', pg_source)
        self.assertIn('Test-Workbook2', pg_source)

    # test the join workgroup link works
    def test_join_workgroup_link(self):
        self.driver.find_element_by_id('join-workgroup').click()
        time.sleep(1)
        pg_source = self.driver.page_source
        self.assertIn('Request to Join a Workgroup', pg_source)


if __name__ == '__main__':
    unittest.main()
