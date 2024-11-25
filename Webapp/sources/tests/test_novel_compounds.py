from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from pony.orm import db_session, select
from sources import app, db
import time


class NovelCompoundLiveTestCase(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """ Need to test the process of adding a novel compound to the reaction database via the sketcher"""

    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        # close the browser window
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        setup_selenium(self)
        login(self, "BB_Test", "BB_login")
        select_workgroup(self)
        select_workbook(self, 0)
        make_new_reaction(self)

    def test_input_fields_auto_fill(self):
        # fill in sketcher, submit, and then assert input fields appear with data
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(5)
        # do the name check because the api call might just be slow
        check_name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        if not check_name:
            self.driver.find_element_by_id('js-new-compound-name').send_keys('pyrazino[2,3-g]quinoxaline-5,10-dione')
        name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        mw = self.driver.find_element_by_id('js-new-compound-mw').get_attribute('value')
        smiles = self.driver.find_element_by_id('js-new-compound-smiles').get_attribute('value')
        self.assertEqual('pyrazino[2,3-g]quinoxaline-5,10-dione', name)
        self.assertEqual('212.03', mw)
        self.assertEqual('O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2', smiles)

    def test_input_fields_accept_input(self):
        # fill in sketcher, submit, and then assert input fields can be filled in
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        self.driver.find_element_by_id('js-new-compound-name').send_keys('TestPhos')
        self.driver.find_element_by_id('js-new-compound-hazards').send_keys('H301')
        self.driver.find_element_by_id('js-new-compound-mw').send_keys('600')
        self.driver.find_element_by_id('js-new-compound-cas').send_keys('123-45-6')
        self.driver.find_element_by_id('js-new-compound-density').send_keys('0.8')
        self.driver.find_element_by_id('js-new-compound-concentration').send_keys('4.0')
        a = self.driver.find_element_by_id('js-new-compound-smiles').is_displayed()
        self.assertFalse(a)

    def test_hazards_valid(self):
        # tests the input fields reject data based on validators
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        check_name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        if not check_name:
            self.driver.find_element_by_id('js-new-compound-name').send_keys('pyrazino[2,3-g]quinoxaline-5,10-dione')
        self.driver.find_element_by_id('js-new-compound-hazards').send_keys('H301, H309')
        self.driver.find_element_by_id('js-new-compound-submit').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        expected_feedback = 'Hazard code "H301, H309" is invalid. Must be valid hazard code and formatted correctly. e.g., H200-H301.'
        self.assertEqual(expected_feedback, alert.text)

    def test_cas_valid(self):
        # tests the input fields reject data based on validators
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        check_name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        if not check_name:
            self.driver.find_element_by_id('js-new-compound-name').send_keys('pyrazino[2,3-g]quinoxaline-5,10-dione')
        self.driver.find_element_by_id('js-new-compound-cas').send_keys('123-45-68-9')
        self.driver.find_element_by_id('js-new-compound-submit').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        expected_feedback = 'CAS invalid.'
        self.assertEqual(expected_feedback, alert.text)

    def test_numerical_inputs_valid(self):
        # tests the input fields reject data based on validators
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        check_name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        if not check_name:
            self.driver.find_element_by_id('js-new-compound-name').send_keys('pyrazino[2,3-g]quinoxaline-5,10-dione')
        self.driver.find_element_by_id('js-new-compound-density').send_keys('not a number')
        self.driver.find_element_by_id('js-new-compound-submit').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        expected_feedback = 'Molecular weight, density, and concentration must be empty or a positive number'
        self.assertEqual(expected_feedback, alert.text)

    def test_unique_name(self):
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        self.driver.find_element_by_id('js-new-compound-name').send_keys('Starting material')
        self.driver.find_element_by_id('js-new-compound-submit').click()
        time.sleep(1)
        alert = self.driver.switch_to.alert
        expected_feedback = 'A compound with this name is already in the database'
        self.assertEqual(expected_feedback, alert.text)

    @db_session
    def test_reactant_added_to_db(self):
        # tests compound is added to the database with the data from the input fields
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        check_name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        if not check_name:
            self.driver.find_element_by_id('js-new-compound-name').send_keys('pyrazino[2,3-g]quinoxaline-5,10-dione')
        self.driver.find_element_by_id('js-new-compound-cas').send_keys('123-45-6')
        self.driver.find_element_by_id('js-new-compound-submit').click()
        time.sleep(1)
        c_name, c_cas = select((c.name, c.cas) for c in db.NovelCompound if c.smiles == 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2' and c.cas == '123-45-6').first()
        self.assertEqual('pyrazino[2,3-g]quinoxaline-5,10-dione', c_name)
        self.assertEqual('123-45-6', c_cas)

    @db_session
    def test_product_added_to_db(self):
        # tests compound is added to the database with the data from the input fields
        rxn_smiles = 'C1=CC=CC=C1>>O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        check_name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        if not check_name:
            self.driver.find_element_by_id('js-new-compound-name').send_keys('pyrazino[2,3-g]quinoxaline-5,10-dione')
        self.driver.find_element_by_id('js-new-compound-cas').send_keys('123-45-6')
        self.driver.find_element_by_id('js-new-compound-submit').click()
        time.sleep(1)
        c_name, c_cas = select((c.name, c.cas) for c in db.NovelCompound if c.smiles == 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2' and c.cas == '123-45-6').first()
        self.assertEqual('pyrazino[2,3-g]quinoxaline-5,10-dione', c_name)
        self.assertEqual('123-45-6', c_cas)

    def test_compound_available_for_wb(self):
        # tests the novel compound added is available for users in that workbook
        rxn_smiles = 'O=C(NC1=C2C=CC=CC2=CC2=C1C=CC=C2)C1=CC=CC=C1>>CCNC(=O)C1=CC=CC=C1.CC'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        novel_compound_name = self.driver.find_element_by_id('js-reactant1').get_attribute('value')
        self.assertIn('Starting material', novel_compound_name)

    def test_compound_not_available_for_other_wbs(self):
        # tests the novel compound added is not available for users in other workbooks
        rxn_smiles = 'CC.OC(=O)C1=CC=CC=C1>>C1CC1C1=CN(N=C1)N1CCNCC1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(4)
        # compound not found so input fields should appear
        mw = self.driver.find_element_by_id('js-new-compound-mw').get_attribute('value')
        self.assertIn('192.14', mw)

    def test_compound_information_generated(self):
        # test the novel compound data is generated for saving into the database - inchi
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        check_name = self.driver.find_element_by_id('js-new-compound-name').get_attribute('value')
        if not check_name:
            self.driver.find_element_by_id('js-new-compound-name').send_keys('pyrazino[2,3-g]quinoxaline-5,10-dione')
        self.driver.find_element_by_id('js-new-compound-submit').click()
        time.sleep(2)
        with db_session:
            c_inchi = select(c.InChI for c in db.NovelCompound if c.smiles == 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2').first()
        self.assertEqual('InChI=1S/C10H4N4O2/c15-9-5-6(12-2-1-11-5)10(16)8-7(9)13-3-4-14-8/h1-4H', c_inchi)

    def test_add_new_reagent_by_table(self):
        """ Test adding a new reagent by table and asserting name appears in the table after"""
        demo_reaction(self)
        self.driver.find_element_by_id("js-add-new-reagent-by-table").click()
        time.sleep(1)
        self.driver.find_element_by_id("js-input-reagent-name").send_keys('unittest_reagent1')
        self.driver.find_element_by_id('js-input-reagent-mw').send_keys(300)
        self.driver.find_element_by_id("js-input-reagent").click()
        time.sleep(1)
        alert = self.driver.switch_to.alert.accept()
        time.sleep(2)
        reagent_name = self.driver.find_element_by_id('js-reagent1').get_attribute('value')
        self.assertEqual('unittest_reagent1', reagent_name)


class NovelCompoundsTestCase(flask_testing.TestCase):
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

    def add_novel_compound(self, name, workbook, workgroup,
                           density='', concentration='', mol_weight='', cas='', smiles='', hphrase=''):
        """Posts data to the _novel_compound routes"""
        return self.client.post(
            '/_novel_compound',
            data=dict(name=name, workbook=workbook, workgroup=workgroup, density=density, concentration=concentration,
                      molWeight=mol_weight, cas=cas, smiles=smiles, hPhrase=hphrase, component="reagent")
        )

    def test_add_novel_compound(self):
        """Tests a novel compound is added"""
        response = self.add_novel_compound(name='test_compound1', workbook='Test-Workbook5', workgroup='Test-Workgroup-3')
        self.assertIn(b"Compound added to the database", response.data)

    def test_duplicate_cas(self):
        """Tests a duplicate cas is detected and compound not added to database"""
        response = self.add_novel_compound(name='test_compound1', workbook='Test-Workbook5', workgroup='Test-Workgroup-3', cas='67-56-1')
        self.assertIn(b'CAS already in database. Please add this compound to the reaction table by searching for the '
                      b'CAS in the reagent box', response.data)

    def test_invalid_smiles(self):
        """ Tests an invalid smiles is detected and compound not added to database"""
        response = self.add_novel_compound(name='test_compound1', workbook='Test-Workbook5', workgroup='Test-Workgroup-3', smiles='not a smiles')
        self.assertIn(b'Invalid smiles', response.data)

    def test_non_numeric_values(self):
        """Tests an entry which is expected to be a positive number"""
        response = self.add_novel_compound(name='test_compound1', workbook='Test-Workbook5', workgroup='Test-Workgroup-3', density='-100')
        self.assertIn(b'Molecular weight, density, and concentration must be empty or a positive number', response.data)

    def test_non_unique_name(self):
        """Tests non-unique name"""
        response = self.add_novel_compound(name='methanol', workbook='Test-Workbook5',
                                           workgroup='Test-Workgroup-3')
        self.assertIn(b'A compound with this name is already in the database', response.data)

    def test_hazard_codes_format(self):
        """Tests hazard code format"""
        response = self.add_novel_compound(name='test_compound1', workbook='Test-Workbook5',
                                           workgroup='Test-Workgroup-3', hphrase='H301, and H320')
        self.assertIn(b'Hazard code \\"H301, and H320\\" is invalid. Must be valid hazard code and formatted correctly. e.g., H200-H301.', response.data)

    def test_add_novel_compound_all_fields(self):
        response = self.add_novel_compound(name='test_compound1', workbook='Test-Workbook5', workgroup='Test-Workgroup-3',
                                           density='1.1', concentration='2.0', mol_weight='300', cas='123-45-6',
                                           smiles='C=CC(=O)N1CCC[C@H](C1)N2C3=C(C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=NC=N3)N',
                                           hphrase='H301-H201-H200-H310')
        self.assertIn(b'Compound added to the database', response.data)

    @db_session
    def test_removes_whitespace(self):
        self.add_novel_compound(name='test       compound1     ', workbook='Test-Workbook5',
                                           workgroup='Test-Workgroup-3')
        name = select(x.name for x in db.NovelCompound if x.name == 'test compound1')
