from auxiliary_for_tests import *
from unittest import mock
from sources import app, db
from pony.orm import db_session, select
from selenium.webdriver.support.select import Select
import flask_testing
from database_setup import test_database_create
import time


@mock.patch('sources.auxiliary.current_user')
@mock.patch('sources.create_workgroup.routes.current_user')
class CreateWorkgroupTestCase(flask_testing.TestCase):
    """
    The tests ensure a workgroup can be created if certain criteria are met and that it is made correctly.
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

    def create_workgroup(self, institution, workgroup, info="I run a big lab"):
        """Posts data to the _save_reaction routes"""
        return self.client.post(
            '/create_workgroup',
            data=dict(institution=institution, workgroup=workgroup, info=info)
        )

    @db_session
    def test_unique_workgroup(self, mock_user, mock_user1):
        """Test when a workgroup with a unique name it is made"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workgroup("Test University", "New workgroup")
        wg_name = select(wg.name for wg in db.WorkGroup_request if wg.name == 'New workgroup'
                         and "PI@test.com" in wg.principal_investigator.user.email).first()
        self.assertEqual("New workgroup", wg_name)

    @db_session
    def test_non_unique_workgroup(self, mock_user, mock_user1):
        """Test when a workgroup with a non-unique name within an institution, requests are not duplicated"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workgroup("Test University", "New workgroup")
        wg_name = select(wg.name for wg in db.WorkGroup_request if wg.name == 'New workgroup'
                         and "PI@test.com" in wg.principal_investigator.user.email).first()
        self.assertEqual("New workgroup", wg_name)
        # Repeat same test workgroup creation and it should not be made.
        response = self.create_workgroup("Test University", "New workgroup")
        wgs = select(wg.id for wg in db.WorkGroup_request if wg.name == 'New workgroup'
                     and "PI@test.com" in wg.principal_investigator.user.email)[:]
        self.assertEqual(1, len(wgs))

    """ uncomment out when institutions return """
    # @db_session
    # def test_non_unique_different_workgroup(self, mock_user, mock_user1):
    #     """Test a workgroup with the same name as a workgroup at a different institution can be made"""
    #     wgs = select(wg.id for wg in db.WorkGroup_request if wg.name == 'New workgroup'
    #                  and "PI@test.com" in wg.principal_investigator.user.email)[:]
    #     self.assertEqual(0, len(wgs))
    #     mock_user.email = 'PI@test.com'
    #     mock_user1.email = 'PI@test.com'
    #     response = self.create_workgroup("Test University", "New workgroup")
    #     wgs = select(wg.id for wg in db.WorkGroup_request if wg.name == 'New workgroup'
    #                  and "PI@test.com" in wg.principal_investigator.user.email)[:]
    #     self.assertEqual(1, len(wgs))
    #     response2 = self.create_workgroup("University of Testing", "New workgroup")
    #     wgs = select(wg.id for wg in db.WorkGroup_request if wg.name == 'New workgroup'
    #                  and "PI@test.com" in wg.principal_investigator.user.email)[:]
    #     self.assertEqual(2, len(wgs))

    """ uncomment out when institutions return """
    # @db_session
    # def test_invalid_institution_fail(self, mock_user, mock_user1):
    #     """Test an invalid institution does not validate and does not make a new workgroup"""
    #     mock_user.email = 'PI@test.com'
    #     mock_user1.email = 'PI@test.com'
    #     response = self.create_workgroup("Made up University", "New workgroup")
    #     wgs = select(wg.id for wg in db.WorkGroup_request if wg.name == 'New workgroup'
    #                  and "PI@test.com" in wg.principal_investigator.user.email)[:]
    #     self.assertEqual(0, len(wgs))

    @db_session
    def test_current_user_made_PI(self, mock_user, mock_user1):
        """This function tests the currently logged in user who makes the new group becomes a PI of that group"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        response = self.create_workgroup("Test University", "New workgroup2")
        pi_email = select(wg.principal_investigator.user.email for wg in db.WorkGroup
                          if 'PI@test.com' in wg.principal_investigator.user.email).first()
        self.assertEqual("PI@test.com", pi_email)

    @db_session
    def test_case_insensitive(self, mock_user, mock_user1):
        """Tests workgroup name is not case sensitive"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        self.create_workgroup("Test University", "New workgroup")
        wgs_before = select(wg.id for wg in db.WorkGroup_request)[:]
        # Repeat same test workgroup creation but with changed case and it should not be made.
        self.create_workgroup("Test University", "nEw WorkGroup")
        wgs_after = select(wg.id for wg in db.WorkGroup_request)[:]
        self.assertEqual(len(wgs_before), len(wgs_after))

    @db_session
    def test_whitespace_removed(self, mock_user, mock_user1):
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        self.create_workgroup("Test University", "New workgroup")
        wgs_before = select(wg.id for wg in db.WorkGroup_request)[:]
        # Repeat same test workgroup creation but with changed whitespace and it should not be made.
        self.create_workgroup("Test University", "      New       workgroup      ")
        wgs_after = select(wg.id for wg in db.WorkGroup_request)[:]
        self.assertEqual(len(wgs_before), len(wgs_after))
        wg_name = select(wg.name for wg in db.WorkGroup_request if wg.name == 'New workgroup'
                     and "PI@test.com" in wg.principal_investigator.user.email).first()
        self.assertEqual('New workgroup', wg_name)

    @db_session
    def test_alphanumeric(self, mock_user, mock_user1):
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        wgs_before = select(wb for wb in db.WorkGroup_request)[:]
        self.create_workgroup("Test University", "new &&&&***@@@_workgroup")
        wgs_after = select(wb for wb in db.WorkGroup_request)[:]
        self.assertEqual(len(wgs_before), len(wgs_after))

    @db_session
    def test_approve_status(self, mock_user, mock_user1):
        """Test the approval status is set to false when workgroup is made"""
        mock_user.email = 'PI@test.com'
        mock_user1.email = 'PI@test.com'
        self.create_workgroup("Test University", "New workgroup")
        wg_status = select(wg.approved for wg in db.WorkGroup if wg.name == 'New workgroup').first()
        self.assertEqual(False, wg_status)


class CreateWorkgroupBrowserTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):

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
        self.driver.find_element_by_id('create-workgroup').click()

    def create_workgroup(self, workgroup, info_text='I am a PI'):
        clear_and_send_keys(self, 'workgroup', workgroup)
        # Fill the textbox for new workgroup name and PI message
        clear_and_send_keys(self, 'info', info_text)
        self.driver.find_element_by_id('submit').click()

    def test_workgroup_creation_approval(self):
        """Here we test clicking from the homepage to the workgroup creation page, enter a unique workgroup name and
        confirm it has been made"""
        # Fill the textbox for new workgroup name and PI message
        self.create_workgroup('New workgroup')
        time.sleep(2)
        # Confirm success message is rendered on page
        self.assertIn("Your Workgroup has been created. It will show as under moderation until it has been approved by "
                      "the site admin.You will receive a notification when your request has been approved.",
                      self.driver.page_source)
        # confirm workgroup is in dropdown
        select_wg = Select(self.driver.find_element_by_id("WG-select"))
        select_wg.select_by_visible_text('New workgroup')
        self.driver.find_element_by_id('go-to-workgroup').click()
        # confirm moderation warning and name are present
        self.assertIn('<h2>New workgroup</h2>\n     <h6 style="color: red">Workgroup Pending Moderation</h6>',
                      self.driver.page_source)
        time.sleep(1)
        # log in as admin
        logout(self)
        login(self, "admin", "admin_login")
        # check notification
        check_notifications(self)
        self.assertIn("A PI has requested a new Workgroup. Go to the &lt;a href=/admin_dashboard id=&#34;admin_"
                      "dashboard&#34;&gt;admin dashboard&lt;/a&gt; to see this request.", self.driver.page_source)
        self.driver.find_element_by_id('TopNavHomeButton').click()
        # approve request
        self.driver.find_element_by_id('admin_dashboard').click()
        self.driver.find_element_by_id("request-approve1").click()
        self.assertIn("This request has been approved", self.driver.page_source)
        # log in as PI again
        logout(self)
        login(self, "PI_Test", "PI_login")
        # check notification
        check_notifications(self)
        self.assertIn("has been approved", self.driver.page_source)
        # check access to new workgroup
        # Confirm workgroup has been added to database
        with db_session:
            a = select(wg.name for wg in db.WorkGroup if wg.name == 'New workgroup'
                       and "PI@test.com" in wg.principal_investigator.user.email).first()
        self.assertIn("New workgroup", a)
        # confirm workgroup is in dropdown
        self.driver.find_element_by_id('TopNavHomeButton').click()
        select_wg = Select(self.driver.find_element_by_id("WG-select"))
        select_wg.select_by_visible_text('New workgroup')
        self.driver.find_element_by_id('go-to-workgroup').click()
        # confirm moderation warning has been removed
        self.assertNotIn('<h6 style="color: red">Workgroup Pending Moderation</h6>',
                      self.driver.page_source)

    def test_workgroup_creation_denial(self):
        """Here we test clicking from the homepage to the workgroup creation page, enter a unique workgroup name and
        confirm it has been made"""
        # Fill the textbox for new workgroup name and PI message
        self.create_workgroup('New workgroup 2')
        time.sleep(1)
        # Confirm success message is rendered on page
        self.assertIn(
            "Your Workgroup has been created. It will show as under moderation until it has been approved by the site admin."
            "You will receive a notification when your request has been approved.",
            self.driver.page_source)
        # log in as admin
        logout(self)
        login(self, "admin", "admin_login")
        # check notification
        check_notifications(self)
        self.assertIn("A PI has requested a new Workgroup. Go to the &lt;a href=/admin_dashboard id=&#34;admin_"
                      "dashboard&#34;&gt;admin dashboard&lt;/a&gt; to see this request.", self.driver.page_source)
        self.driver.find_element_by_id('TopNavHomeButton').click()
        # approve request
        self.driver.find_element_by_id('admin_dashboard').click()
        self.driver.find_element_by_id("request-deny1").click()
        alert = self.driver.switch_to.alert
        self.assertIn("Press OK if you are sure you mean to delete the workgroup:New workgroup 2", alert.text)
        alert.accept()
        self.assertIn("This request has been denied", self.driver.page_source)
        # log in as PI again
        logout(self)
        login(self, "PI_Test", "PI_login")
        time.sleep(1)
        # check notification
        check_notifications(self)
        self.assertIn("has been denied", self.driver.page_source)
        # check access to new workgroup
        # Confirm workgroup has been added to database
        with db_session:
            a = select(wg.name for wg in db.WorkGroup if wg.name == 'New workgroup 2'
                       and "PI@test.com" in wg.principal_investigator.user.email).first()
        self.assertIsNone(a)

    def test_workgroup_failed_creation(self):
        """Testing a fail message is sent when a non-unique workgroup is made"""
        with db_session:
            # Add a workgroup to the database directly and then check we can't make another with the same name
            institution = select(i for i in db.Institution if i.name == 'University of Testing').first()
            pi = select(p for p in db.Person if p.user.email == 'PI@test.com').first()
            wg = db.WorkGroup(name="New workgroup 10", institution=institution, principal_investigator=pi)
        # Select an institution on the dropdown
        # institution = Select(self.driver.find_element_by_id("institution"))
        # institution.select_by_visible_text("University of Testing")
        # Fill the textbox for new workgroup name
        self.create_workgroup('New workgroup 10')
        # Confirm success message is rendered on page
        self.assertIn("A Workgroup of this name already exists", self.driver.page_source)
        # Confirm workgroup has been added to database
        with db_session:
            wgs = select(wg.id for wg in db.WorkGroup if wg.name == 'New workgroup 10'
                         and "PI@test.com" in wg.principal_investigator.user.email)[:]
        self.assertEqual(1, len(wgs))

    def test_only_one_unapproved_workgroup(self):
        """Testing a fail message is sent when the creator already has 1 unapproved workgroup"""
        # Fill the textbox for new workgroup name and PI message
        self.create_workgroup('New workgroup')
        time.sleep(2)
        # Confirm success message is rendered on page
        self.assertIn(
            "Your Workgroup has been created. It will show as under moderation until it has been approved by the site admin."
            "You will receive a notification when your request has been approved.",
            self.driver.page_source)
        # navigate back to create workgroup and try to create another workgroup
        self.driver.find_element_by_id('create-workgroup').click()
        time.sleep(1)
        self.create_workgroup('New new workgroup')
        # Confirm success message is rendered on page
        self.assertIn("You cannot create another workgroup until your first workgroup has been approved", self.driver.page_source)
        # Confirm workgroup has not been added to database
        with db_session:
            wgs = select(wg.id for wg in db.WorkGroup if wg.name == 'New new workgroup'
                         and "PI@test.com" in wg.principal_investigator.user.email).first()
        self.assertEqual(None, wgs)
    # def test_institutions_displayed(self):
    #     """This function checks both the institutions added to the database are present in the dropdown"""
    #     institution = Select(self.driver.find_element_by_id("institution"))
    #     name_ls = [op.text for op in institution.options]
    #     self.assertIn("University of Testing", name_ls)
    #     self.assertIn("Test University", name_ls)
