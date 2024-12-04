from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from sources import app, db
from pony.orm import select, db_session
import time


class TestNovelSolvents(flask_testing.LiveServerTestCase):
	def create_app(self):
		app.config.from_object('config.TestConfig')
		# restore_solvent_table_db()
		with db_session:
			select(x for x in db.Solvent if x.name == 'ionic liquid').delete()
		test_database_create(db)
		return app

	def tearDown(self):
		self.driver.quit()
		# restore_solvent_table_db()
		with db_session:
			select(x for x in db.Solvent if x.name == 'ionic liquid').delete()
		restore_db()

	def setUp(self):
		"""We load a test instance of the app, clear the db, register a new user then login"""
		setup_selenium(self)
		login(self)
		select_workgroup(self)
		select_workbook(self)
		make_new_reaction(self)
		demo_reaction(self)

	@db_session
	def newly_added_solvent_assertions(self, solvent_name='ionic liquid', table_number='3'):
		# compound should now be in db
		db_check = select(
			x for x in db.NovelCompound if x.solvent.name.lower() == solvent_name.lower()
			and x.name.lower() == solvent_name.lower()).first()
		self.assertNotEqual(None, db_check)
		# compound should be in table
		table_name_check = self.driver.find_element_by_id("js-solvent1").get_attribute('value')
		self.assertEqual(solvent_name, table_name_check)
		# table number should be 3
		table_number_check = self.driver.find_element_by_id("js-solvent-table-number1").get_attribute('value')
		self.assertEqual(table_number, table_number_check)
		# compound should be in datalist
		datalist = self.driver.find_element_by_id("js-solvent-datalist1")
		# Find all of the option elements within the datalist
		options = datalist.find_elements_by_tag_name("option")
		# Iterate through the options and check their values
		option_found = False
		for option in options:
			if option.get_attribute("value") == "ionic liquid":
				option_found = True
		self.assertEqual(True, option_found)

	# check datalist membership
	@db_session
	def add_novel_solvent_to_db(self):
		workbook = select(x for x in db.WorkBook if x.name == 'Test-Workbook2').first()
		nc = db.NovelCompound(name='ionic liquid', cas='123-45-6', hphrase='Unknown', workbook=workbook)
		db.Solvent(name='ionic liquid', flag=5, hazard='Unknown', compound=None, novel_compound=nc)

	def test_add_novel_solvent_button(self):
		"""Tests adding a novel solvent by the button and checks presence in db"""
		self.driver.find_element_by_id("js-add-new-solvent-by-table").click()
		time.sleep(1)
		self.driver.find_element_by_id("js-input-solvent-name").send_keys('ionic liquid')
		self.driver.find_element_by_id("js-input-solvent").click()
		time.sleep(0.5)
		self.driver.switch_to.alert.accept()
		time.sleep(0.5)
		self.driver.switch_to.alert.accept()
		time.sleep(1)
		self.newly_added_solvent_assertions()

	def test_fields_add_novel_solvent_cas_entered(self):
		"""Tests the submission fields for a novel solvent appear when an unknown cas is entered in the solvent
		search box, & correct table number"""
		self.driver.find_element_by_class_name('js-add-solvent').click()
		self.driver.find_element_by_id("js-solvent1").send_keys("123-45-6")
		time.sleep(1)
		self.driver.switch_to.alert.accept()
		time.sleep(0.5)
		self.driver.find_element_by_id("js-input-solvent-name").send_keys('ionic liquid')
		# cas should automatically be inputted
		cas_check = self.driver.find_element_by_id("js-input-solvent-cas").get_attribute('value')
		self.assertEqual('123-45-6', cas_check)
		self.driver.find_element_by_id("js-input-solvent").click()
		time.sleep(1)
		self.driver.switch_to.alert.accept()
		time.sleep(1)
		self.driver.switch_to.alert.accept()
		time.sleep(1)
		self.newly_added_solvent_assertions()

	def test_colouring_of_novel_solvents(self):
		"""Tests the colouring changes back when switching between novel/non-novel solvents"""
		self.add_novel_solvent_to_db()
		self.driver.refresh()
		time.sleep(2)
		self.driver.find_element_by_id("action-button-submit").click()
		time.sleep(5)
		self.driver.find_element_by_class_name('js-add-solvent').click()
		self.driver.find_element_by_id("js-solvent1").send_keys("ben")
		move_down_one_option_action_chains(self)
		time.sleep(1)
		value = self.driver.find_element_by_id('js-solvent1').get_attribute("value")
		self.assertEqual(value, "Benzene")
		rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
		self.assertEqual(rgb, "rgba(139, 0, 0, 1)")
		clear_and_send_keys(self, "js-solvent1", "ionic liq")
		move_down_one_option_action_chains(self)
		time.sleep(1)
		value = self.driver.find_element_by_id('js-solvent1').get_attribute("value")
		self.assertEqual(value, "ionic liquid")
		rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
		self.assertEqual(rgb, "rgba(255, 255, 255, 1)")

	def test_novel_solvent_appears_in_workbooks_dropdown(self):
		"""Tests the novel solvent appears in solvent dropdown"""
		self.add_novel_solvent_to_db()
		self.driver.refresh()
		time.sleep(2)
		self.driver.find_element_by_id("action-button-submit").click()
		time.sleep(5)
		self.driver.find_element_by_class_name('js-add-solvent').click()
		# Locate the datalist element
		datalist = self.driver.find_element_by_id("js-solvent-datalist1")
		# Find all of the option elements within the datalist
		options = datalist.find_elements_by_tag_name("option")
		# Iterate through the options and check their values
		option_found = False
		for option in options:
			if option.get_attribute("value") == "ionic liquid":
				option_found = True
		self.assertEqual(True, option_found)
		self.driver.find_element_by_id("js-solvent1").send_keys("ionic liq")
		move_down_one_option_action_chains(self)
		value = self.driver.find_element_by_id('js-solvent1').get_attribute("value")
		self.assertEqual(value, "ionic liquid")
		rgb = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
		self.assertEqual(rgb, "rgba(255, 255, 255, 1)")
		class_name = self.driver.find_element_by_id("js-solvent1").get_attribute("class")
		self.assertEqual("non-chem21 remove-highlight-filled-cell", class_name)

	def test_novel_solvent_doesnt_appear_outside_workbooks_dropdown(self):
		"""Tests the novel solvent doesnt appear for users outside the workbook"""
		self.add_novel_solvent_to_db()
		self.driver.find_element_by_id("TopNavHomeButton").click()
		time.sleep(1)
		select_workgroup(self)
		select_workbook(self, 0)
		make_new_reaction(self)
		demo_reaction(self)
		self.driver.find_element_by_class_name('js-add-solvent').click()
		# Locate the datalist element
		datalist = self.driver.find_element_by_id("js-solvent-datalist1")
		# Find all of the option elements within the datalist
		options = datalist.find_elements_by_tag_name("option")
		# Iterate through the options and check their values
		option_found = False
		for option in options:
			if option.get_attribute("value") == "ionic liquid":
				option_found = True
		self.assertEqual(False, option_found)
