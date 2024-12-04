from auxiliary_for_tests import *
from unittest import mock, main
from sources import app, db
from pony.orm import db_session, select
import flask_testing
from database_setup import test_database_create
import time


@mock.patch('sources.auxiliary.current_user')
@mock.patch('sources.create_workbook.routes.current_user')
class CreateWorkbookTestCase(flask_testing.TestCase):
    """
    The tests ensure a workbook can be created if certain criteria are met and that it is made correctly.
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
        # Don't want to drop compound or solvent tables.
        restore_db()

    def create_workbook(self, workbook, workbook_abbreviation, workgroup="Test-Workgroup"):
        """Posts data to the create_workbook routes"""
        return self.client.post(
            '/create_workbook/' + workgroup,
            data=dict(workgroup=workgroup, abbreviation=workbook_abbreviation, workbook=workbook)
            )

    @db_session
    def test_unique_workbook(self, mock_user, mock_user1):
        """Test when a workbook with a unique name it is made"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workbook("New workbook", 'NW1')
        wb_name = select(wb.name for wb in db.WorkBook if wb.name == 'New workbook'
                         and "PI@test.com" in wb.users.user.email).first()
        self.assertIn("New workbook", wb_name)

    @db_session
    def test_non_unique_workbook(self, mock_user, mock_user1):
        """Test when a workbook with a non-unique name within an workgroup, it is not made"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workbook("New workbook 2", 'NW1')
        wb_name = select(wb.name for wb in db.WorkBook if wb.name == 'New workbook 2'
                         and "PI@test.com" in wb.users.user.email).first()
        self.assertIn("New workbook 2", wb_name)
        mock_user.email = 'PI@test.com'
        response = self.create_workbook("New workbook 2", 'NW2')
        wbs = select(wb.id for wb in db.WorkBook if wb.name == 'New workbook 2'
                         and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(1, len(wbs))

    @db_session
    def test_non_unique_workbook_abbreviation(self, mock_user, mock_user1):
        """Test when a workbook with a non-unique name within an workgroup, it is not made"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workbook("New workbook 2", 'NW1')
        wb_name = select(wb.name for wb in db.WorkBook if wb.name == 'New workbook 2'
                         and "PI@test.com" in wb.users.user.email).first()
        self.assertIn("New workbook 2", wb_name)
        mock_user.email = 'PI@test.com'
        response = self.create_workbook("New workbook 3", 'NW2')
        wbs = select(wb.id for wb in db.WorkBook if wb.name == 'New workbook 2'
                         and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(1, len(wbs))

    @db_session
    def test_non_unique_workbook_different_workgroups(self, mock_user, mock_user1):
        """Test a workbook with the same name as a workbook in a different workgroup can be made"""
        wbs = select(wb.id for wb in db.WorkBook if wb.name == 'New workbook'
                         and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(0, len(wbs))
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workbook("New workbook", 'NW1')
        wbs = select(wb.id for wb in db.WorkBook if wb.name == 'New workbook'
                         and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(1, len(wbs))
        response = self.create_workbook("New workbook", 'NW2', workgroup='Test-Workgroup-2')
        wbs = select(wb.id for wb in db.WorkBook if wb.name == 'New workbook'
                         and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(2, len(wbs))

    @db_session
    def test_current_user_made_PI(self, mock_user, mock_user1):
        """This function tests the currently logged in user who makes the new workbook becomes a user of that workbook"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workbook("New workbook", 'NW1')
        pi_email = select(wb.users.user.email for wb in db.WorkBook if wb.name == 'New workbook'
                         and "PI@test.com" in wb.users.user.email).first()
        self.assertEqual('PI@test.com', pi_email)

    @db_session
    def test_current_user_made_SR(self, mock_user, mock_user1):
        mock_user.email = 'SR@test.com'
        mock_user1.email = 'SR@test.com'
        response = self.create_workbook("New workbook", 'NW1')
        sr_email = select(wb.users.user.email for wb in db.WorkBook if wb.name == 'New workbook'
                         and "SR@test.com" in wb.users.user.email).first()
        self.assertEqual('SR@test.com', sr_email)

    @db_session
    def test_case_insensitive(self, mock_user, mock_user1):
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workbook("new workbook", 'NW1')
        wbs = select(wb.id for wb in db.WorkBook if wb.name == 'new workbook'
                         and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(1, len(wbs))
        response = self.create_workbook("NeW WorKBook", 'NW1')
        wbs = select(wb.id for wb in db.WorkBook if wb.name == 'new workbook'
                         and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(1, len(wbs))

    @db_session
    def test_whitespace_removed(self, mock_user, mock_user1):
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workbook("    new        workbook     ", 'NW1')
        wb_name = select(wb.name for wb in db.WorkBook if wb.name == 'new workbook'
                         and 'PI@test.com' in wb.users.user.email).first()
        self.assertEqual('new workbook', wb_name)

    @db_session
    def test_alphanumeric(self, mock_user, mock_user1):
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        wbs_before = select(wb for wb in db.WorkBook)[:]
        response = self.create_workbook("new &&&&***@@@_workbook", 'NW1')
        wbs_after = select(wb for wb in db.WorkBook)[:]
        self.assertEqual(len(wbs_before), len(wbs_after))


class CreateWorkbookBrowserTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):

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
        select_workgroup(self)
        self.driver.find_element_by_id("manage-workbook").click()
        self.driver.find_element_by_id("create-workbook-button").click()

    def test_workbook_creation(self):
        """Here we test clicking from the homepage to the workbook creation page, enter a unique workbook name and
        confirm it has been made"""
        # Fill the textbox for new workbook name
        clear_and_send_keys(self, "workbook", "New workbook")
        self.driver.find_element_by_id('abbreviation').send_keys('NW1')
        self.driver.find_element_by_id("submit").click()
        # Confirm success message is rendered on page
        time.sleep(1)
        self.assertIn("Workbook has been created", self.driver.page_source)
        # Confirm workbook has been added to database
        with db_session:
            wb_name = select(wb.name for wb in db.WorkBook if wb.name == 'New workbook'
                             and "PI@test.com" in wb.users.user.email).first()
        self.assertIn("New workbook", wb_name)

    def test_workbook_failed_creation(self):
        """Testing a fail message is sent when a non-unique workbook is made"""
        with db_session:
            workgroup = select(wg for wg in db.WorkGroup if "Test-Workgroup" == wg.name).first()
            user = select(p for p in db.Person if p.user.email == 'PI@test.com').first()
            wb = db.WorkBook(name="New workbook", abbreviation='NW1', group=workgroup, users=user)
        # Fill the textbox for new workbook name
        clear_and_send_keys(self, "workbook", "New workbook")
        self.driver.find_element_by_id('abbreviation').send_keys('NW1')
        self.driver.find_element_by_id("submit").click()
        # Confirm success message is rendered on page
        self.assertIn("A workbook of this name already exists", self.driver.page_source)
        # Confirm workbook has been added to database
        with db_session:
            wbs = select(wb.id for wb in db.WorkBook if wb.name == 'New workbook'
                             and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(1, len(wbs))

    def test_workbook_failed_creation_duplicate_abbreviation(self):
        """Testing a fail message is sent when a non-unique workbook is made"""
        with db_session:
            workgroup = select(wg for wg in db.WorkGroup if "Test-Workgroup" == wg.name).first()
            user = select(p for p in db.Person if p.user.email == 'PI@test.com').first()
            wb = db.WorkBook(name="New workbook", abbreviation='NW1', group=workgroup, users=user)
        # Fill the textbox for new workbook name
        clear_and_send_keys(self, "workbook", "New workbook2")
        self.driver.find_element_by_id('abbreviation').send_keys('NW1')
        self.driver.find_element_by_id("submit").click()
        # Confirm success message is rendered on page
        self.assertIn("A workbook with this abbreviation already exists", self.driver.page_source)
        # Confirm workbook has been added to database
        with db_session:
            wbs = select(wb.id for wb in db.WorkBook if wb.name == 'New workbook'
                         and "PI@test.com" in wb.users.user.email)[:]
        self.assertEqual(1, len(wbs))


class CreateWorkbookTestSeniorResearcher(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """This test checks a senior researcher is also able to create a new workbook"""
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
        self.driver.find_element_by_id("manage-workbook").click()
        time.sleep(1)
        self.driver.find_element_by_id("create-workbook-button").click()
        time.sleep(1)

    def test_workbook_creation(self):
        """Here we test clicking from the homepage to the workbook creation page, enter a unique workbook name and
        confirm it has been made"""
        # Fill the textbox for new workbook name
        clear_and_send_keys(self, "workbook", "New workbook")
        self.driver.find_element_by_id('abbreviation').send_keys('NW1')
        self.driver.find_element_by_id("submit").click()
        # Confirm success message is rendered on page
        self.assertIn("Workbook has been created", self.driver.page_source)
        # Confirm workbook has been added to database
        with db_session:
            wb_name = select(wb.name for wb in db.WorkBook if wb.name == 'New workbook'
                             and "SR@test.com" in wb.users.user.email).first()
        self.assertIn("New workbook", wb_name)


if __name__ == '__main__':
    main()
