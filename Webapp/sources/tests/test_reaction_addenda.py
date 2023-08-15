from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create


# test all the buttons are rendered and linked
class ReactionsAddendaTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
	# create the app
	def create_app(self):
		app.config.from_object('config.TestConfig')
		test_database_create(db)
		return app

	def tearDown(self):
		self.driver.quit()
		restore_db()

	def setUp(self):
		"""We load a test instance of the app"""
		setup_selenium(self)
		login(self)
		select_workgroup(self)
		select_workbook(self, idx=0)
		make_and_lock_demo_reaction(self)

	def test_note_can_be_added_after_locking(self):
		"""Tests the 'Add Note' button appears by entering a a new note and checking it appears"""
		self.driver.find_element_by_id("reaction-note-button").click()
		time.sleep(1)
		self.driver.find_element_by_id("new-reaction-note-text").send_keys("Hello world")
		self.driver.find_element_by_id("new-reaction-note-submit").click()
		time.sleep(1)
		addendum_text = self.driver.find_element_by_id("addendum-text1").text
		addendum_author = self.driver.find_element_by_id("addendum-author1").text
		self.assertEqual("Hello world", addendum_text)
		self.assertEqual("Pat Inglis", addendum_author)

	def test_note_appears_after_reloading(self):
		self.driver.find_element_by_id("reaction-note-button").click()
		time.sleep(1)
		self.driver.find_element_by_id("new-reaction-note-text").send_keys("Hello world")
		self.driver.find_element_by_id("new-reaction-note-submit").click()
		self.driver.find_element_by_id("TopNavHomeButton").click()
		time.sleep(1)
		select_workgroup(self)
		select_workbook(self, idx=0)
		self.driver.find_element_by_id('new reaction-reload').click()
		time.sleep(8)
		addendum_text = self.driver.find_element_by_id("addendum-text1").text
		addendum_author = self.driver.find_element_by_id("addendum-author1").text
		self.assertEqual("Hello world", addendum_text)
		self.assertEqual("Pat Inglis", addendum_author)

	def test_symbol_on_reaction_list(self):
		self.driver.find_element_by_id("reaction-note-button").click()
		time.sleep(1)
		self.driver.find_element_by_id("new-reaction-note-text").send_keys("Hello world")
		self.driver.find_element_by_id("new-reaction-note-submit").click()
		self.driver.find_element_by_id("TopNavHomeButton").click()
		time.sleep(1)
		select_workgroup(self)
		select_workbook(self, idx=0)
		pg_source = self.driver.page_source
		self.assertIn("Contains Amendments", pg_source)
		self.assertIn("bi bi-chat-dots", pg_source)
		
	def test_other_user_can_only_see_note_but_not_button(self):
		self.driver.find_element_by_id("reaction-note-button").click()
		time.sleep(1)
		self.driver.find_element_by_id("new-reaction-note-text").send_keys("Hello world")
		self.driver.find_element_by_id("new-reaction-note-submit").click()
		logout(self)
		login(self, "SR_Test", "SR_login")
		select_workgroup(self)
		select_workbook(self, idx=0)
		self.driver.find_element_by_id('new reaction-reload').click()
		time.sleep(8)
		# check element not visible
		reaction_note_button = self.driver.find_element_by_id("reaction-note-button").is_displayed()
		self.assertEqual(False, reaction_note_button)
		addendum_text = self.driver.find_element_by_id("addendum-text1").text
		addendum_author = self.driver.find_element_by_id("addendum-author1").text
		self.assertEqual("Hello world", addendum_text)
		self.assertEqual("Pat Inglis", addendum_author)
		






