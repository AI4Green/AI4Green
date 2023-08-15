from auxiliary_for_tests import *
from unittest import mock
from pony.orm import select, db_session
import flask_testing
from database_setup import test_database_create
from datetime import datetime
from sources import app, db, auxiliary
import time


@mock.patch('sources.save_reaction.routes.current_user')
class SaveNewReactionTestCase(flask_testing.TestCase):
    """
    The tests the process of saving a new reaction
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

    def save_new_reaction(self, name, smiles, reaction_id):
        """Posts data to the _save_reaction routes"""
        return self.client.post(
            '/new_reaction',
            data=dict(reactionID=reaction_id, workgroup="Test-Workgroup", workbook="Test-Workbook",
                      reactionSmiles=smiles, reactionName=name)
        )

    @db_session
    def test_reaction_made(self, mock_user):
        """Tests reaction is made"""
        mock_user.email = 'SM@test.com'
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        reaction = select(
            rxn for rxn in db.Reaction if rxn.name == 'Amide coupling' and rxn.reaction_id == 'TW1-001').first()
        self.assertIsNotNone(reaction)

    @db_session
    def test_reaction_duplicate_id(self, mock_user):
        """Tests reaction is not made with a duplicate id"""
        mock_user.email = 'SM@test.com'
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        reaction = select(
            rxn for rxn in db.Reaction if rxn.name == 'Amide coupling' and rxn.reaction_id == 'TW1-001').first()
        self.assertIsNotNone(reaction)
        response = self.save_new_reaction('Amide coupling2', "CC.CC>>CC.CC", 'TW1-001')
        self.assertIn(b'A reaction with this ID already exists. Please refresh the page and try again.', response.data)

    def test_non_unique_name_within_workbook(self, mock_user):
        """Reaction should not save if name is not unique within a workbook"""
        mock_user.email = 'SM@test.com'
        # Name is unique so the reaction should save
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        self.assertIn(b"New reaction made", response.data)
        # Name is no longer unique and should return warning
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-002')
        self.assertIn(b"This reaction name is already used. Please choose another name.", response.data)

    def test_second_reaction_id(self, mock_user):
        """Reaction should save if the second id and name are unique"""
        mock_user.email = 'SM@test.com'
        # Name is unique so the reaction should save
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        self.assertIn(b"New reaction made", response.data)
        # Name is no longer unique and should return warning
        response = self.save_new_reaction('Amide coupling2', "CC.CC>>CC.CC", 'TW1-002')
        self.assertIn(b"New reaction made", response.data)


@mock.patch('sources.save_reaction.routes.current_user')
class AutoSaveReactionTestCase(flask_testing.TestCase):
    """
    The tests the process of saving a new reaction
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

    def autosave_reaction(self, name, reactants, products, smiles, reaction_id, description='', solvent='',
                          reagents='',
                          complete="not complete", unreacted_product_mass="5", real_product_mass="6"):
        """Posts data to the _save_reaction routes"""
        return self.client.post(
            '/_autosave',
            data=dict(reactionID=reaction_id, workgroup="Test-Workgroup", workbook="Test-Workbook",
                      reactionSmiles=smiles, reactionName=name, reactionDescription=description,
                      reactantPrimaryKeys=reactants, productPrimaryKeys=products, reagentPrimaryKeys=reagents,
                      solventPrimaryKeys=solvent, reactantMasses=['100', '200'], reactantMassesRaw=['100', '200'],
                      reactantAmounts=['3', '5'], reactantAmountsRaw=["3", "5"], reactantEquivalents=['1', '2'],
                      reactantConcentrations=["", ""], reactantVolumes=["", ""], reactantVolumesRaw=["", ""],
                      reactantPhysicalForms=['Dense solid', 'Dense solid'], reactantDensities=["", ""],
                      reactantNames='reactant1;reactant2;', reactantPhysicalFormsText='Solid;Solid;',
                      reactantHazards='H301;H200;', reactantMolecularWeights="100;200;",

                      reagentSmiles='smiles1;smiles2;', reagentPhysicalFormsText='Solid;Solid;',
                      reagentEquivalents=['1', '2'], reagentPhysicalForms=['Dense solid', 'Dense solid'],
                      reagentMolecularWeights=["100", "200"], reagentMasses=["100", "200"],
                      reagentMassesRaw=["100", "200"], reagentVolumes=["", ""], reagentVolumesRaw=["", ""],
                      reagentHazards=["", ""], reagentNames=["", ""],
                      reagentAmounts=["3", "6"], reagentAmountsRaw=["3", "6"],
                      reagentDensities=["", ""], reagentConcentrations=["", ""],

                      solventPhysicalFormsText="Liquid;Liquid;",
                      solventVolumes=['1', '2'], solventNames=["", ""], solventHazards=["", ""],
                      solventConcentrations=["", ""], solventPhysicalForms=['Non-volatile liquid (b.p. > 130 °C',
                                                                            'Volatile liquid (70 °C ≤ b.p. ≤ 130 °C)'],

                      limitingReactantTableNumber='1', mainProductTableNumber='3',
                      productNames='Product1;Product2;', productPhysicalFormsText="Solid;Solid;",
                      productPhysicalForms=['Dense solid', 'Dense solid'], productHazards="H300;H301;",
                      productMolecularWeights="300;400;",

                      amountUnits='mmol', massUnits='mg', volumeUnits='mL', solventVolumeUnits='mL',
                      realProductMass=real_product_mass, productAmounts=["", ""], productAmountsRaw=["", ""],
                      productMasses=["", ""], productMassesRaw=["", ""],
                      unreactedReactantMass=unreacted_product_mass, reactionTemperature=100,
                      elementSustainability='500+ years', batchFlow='Batch', productAmountUnits="mmol",
                      productMassUnits="mg",
                      isolationMethod='Filtration', catalystUsed='No catalyst or reagents', catalystRecovered='-select',
                      selectedRadioButtons=['#HPLC', '#pyrophorics'],
                      customProtocol1='', customProtocol2='', otherHazardTextArea='', researcher="", supervisor="",
                      complete=complete, summary_to_print='no', massEfficiency='90', toExport='no')
        )

    def save_new_reaction(self, name, smiles, reaction_id):
        """Posts data to the _save_reaction routes"""
        return self.client.post(
            '/new_reaction',
            data=dict(reactionID=reaction_id, workgroup="Test-Workgroup", workbook="Test-Workbook",
                      reactionSmiles=smiles, reactionName=name)
        )

    def compound_ids_to_smiles(self, ids):
        expected_ids = ids[:-1].split(';')
        expected_smiles = auxiliary.get_smiles(expected_ids)
        return expected_smiles

    @db_session
    def test_successful_reaction_all_fields(self, mock_user):
        """Test a successful reaction with all of the fields filled"""
        mock_user.email = 'SM@test.com'
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        reaction = select(
            rxn for rxn in db.Reaction if rxn.name == 'Amide coupling' and rxn.reaction_id == 'TW1-001').first()
        self.assertIsNotNone(reaction)
        reactants = '1;2;'
        products = '3;4;'
        reagents = '5;6;'
        solvents = '7;8;'
        description = 'Amide coupling using heterogeneous catalysis'
        response = self.autosave_reaction('Amide coupling', reactants, products, "CC.CC>>CC.CC", 'TW1-001', reagents=reagents,
                                      solvent=solvents, description=description)
        self.assertIn(b'Reaction Updated!', response.data)

    @db_session
    def test_components_saved_correctly(self, mock_user):
        """Check database contents to confirm reaction components have been saved successfully"""
        mock_user.email = 'SM@test.com'
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        reaction = select(
            rxn for rxn in db.Reaction if rxn.name == 'Amide coupling' and rxn.reaction_id == 'TW1-001').first()
        self.assertIsNotNone(reaction)
        reactants = '1;2;'
        products = '3;4;'
        reagents = '5;6;'
        solvents = '7;8;'
        response = self.autosave_reaction('Amide coupling', reactants, products, "CC.CC>>CC.CC", 'TW1-001', reagents=reagents,
                                      solvent=solvents)
        rxn_components = select((x.reactants, x.products, x.reagents, x.solvent) for x in db.Reaction
                                if x.creator.user.email == 'SM@test.com')[:]
        expected_reactants = self.compound_ids_to_smiles(reactants)
        # expected_reagents = self.compound_ids_to_smiles(reagents)
        # expected_solvents = self.compound_ids_to_smiles(solvents)
        expected_products = self.compound_ids_to_smiles(products)
        # smiles for products/reactants processed differently to reagents and solvents.
        self.assertIn((expected_reactants, expected_products, ['smiles1', 'smiles2', ''], ['7', '8']), rxn_components)


    @db_session
    def test_meta_data_saved_correctly(self, mock_user):
        """Checks database contents to confirm reaction meta data has been saved successfully"""
        mock_user.email = 'SM@test.com'
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        reaction = select(
            rxn for rxn in db.Reaction if rxn.name == 'Amide coupling' and rxn.reaction_id == 'TW1-001').first()
        self.assertIsNotNone(reaction)
        reactants = '1;2;'
        products = '3;4;'
        description = 'Amide coupling using heterogeneous catalysis'
        response = self.autosave_reaction('Amide coupling', reactants, products, "CC.CC>>CC.CC", 'TW1-001', description=description)
        rxn_data = list(select([x.name, x.description, x.workbooks.name, x.creator.user.email]
                          for x in db.Reaction if x.creator.user.email == 'SM@test.com')[:])
        rxn_time = select(x.time_of_update for x in db.Reaction
                          if x.creator.user.email == 'SM@test.com').first()
        self.assertEqual([('Amide coupling', 'Amide coupling using heterogeneous catalysis', 'Test-Workbook',
                          'SM@test.com')], rxn_data)
        self.assertIsInstance(rxn_time, datetime)

    @db_session
    def test_update_reaction(self, mock_user):
        mock_user.email = 'SM@test.com'
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        reaction_original = select(
            rxn for rxn in db.Reaction if rxn.name == 'Amide coupling' and rxn.reaction_id == 'TW1-001').first()
        self.assertIsNotNone(reaction_original)
        workbook = select(b for b in db.WorkBook if 'SM@test.com' in b.users.user.email).first()
        # update reaction
        mock_user.email = 'SM@test.com'
        reactants = '1;2;'
        products = '6;7;'
        response = self.autosave_reaction('Amide coupling', reactants, products, "CC.CC>>CC.CC", 'TW1-001')
        self.assertIn(b"Reaction Updated!", response.data)
        # check reaction updated correctly
        rxn_updated, rxn_update_time, rxn_products = select((x, x.time_of_update, x.products) for x in db.Reaction if
                                                       x.creator.user.email == 'SM@test.com' and x.workbooks == workbook).first()
        self.assertEqual(reaction_original, rxn_updated)
        expected_products = products[:-1].split(';')
        expected_product_smiles = auxiliary.get_smiles(expected_products)
        self.assertEqual(rxn_products, expected_product_smiles)
        self.assertIsNotNone(rxn_update_time)

    @db_session
    def test_mark_complete_reaction(self, mock_user):
        mock_user.email = 'SM@test.com'
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        reaction = select(
            rxn for rxn in db.Reaction if rxn.name == 'Amide coupling' and rxn.reaction_id == 'TW1-001').first()
        self.assertIsNotNone(reaction)
        reactants = '1;2;'
        products = '3;4;'
        response = self.autosave_reaction('Amide coupling', reactants, products, "CC.CC>>CC.CC", 'TW1-001', complete="complete")
        # check reaction saved
        workbook = select(b for b in db.WorkBook if 'SM@test.com' in b.users.user.email).first()
        reaction, rxn_id_original = select(
            (x, x.id) for x in db.Reaction if x.creator.user.email == 'SM@test.com' and x.workbooks == workbook and x.complete == "complete").first()
        self.assertIsNotNone(reaction)

    @db_session
    def test_mark_complete_reaction_missing_data(self, mock_user):
        mock_user.email = 'SM@test.com'
        response = self.save_new_reaction('Amide coupling', "CC.CC>>CC.CC", 'TW1-001')
        reaction = select(
            rxn for rxn in db.Reaction if rxn.name == 'Amide coupling' and rxn.reaction_id == 'TW1-001').first()
        self.assertIsNotNone(reaction)
        reactants = '1;2;'
        products = '3;4;'
        response = self.autosave_reaction('Amide coupling', reactants, products, "CC.CC>>CC.CC", 'TW1-001', complete="complete", unreacted_product_mass="", real_product_mass="")
        # check reaction saved
        workbook = select(b for b in db.WorkBook if 'SM@test.com' in b.users.user.email).first()
        try:
            reaction, rxn_id_original = select(
                (x, x.id) for x in db.Reaction if
                x.creator.user.email == 'SM@test.com' and x.workbooks == workbook and x.complete == "complete").first()
        except TypeError:
            reaction = None
        self.assertIsNone(reaction)


class BrowserAutosaveTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """tests the reaction reloads successfully without loss of data even after multiple reloads"""
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
        login(self, "BB_Test", "BB_login")
        select_workgroup(self)
        select_workbook(self, 0)

    def exit_and_reload_reaction(self, reload_id='new reaction-reload'):
        self.driver.find_element_by_id('TopNavHomeButton').click()
        time.sleep(1)
        select_workgroup(self)
        select_workbook(self, 0)
        self.driver.find_element_by_id(reload_id).click()
        time.sleep(6)

    def assert_reactants_data(self):
        mass = self.driver.find_element_by_id("js-reactant-rounded-mass1").get_attribute('value')
        equiv = self.driver.find_element_by_id("js-reactant-equivalent2").get_attribute('value')
        self.assertEqual('100', mass)
        self.assertEqual('2', equiv)

    def assert_smiles_data(self, expected_smiles="OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1"):
        smiles_hidden_elem = self.driver.find_element_by_id('js-reaction-smiles').get_attribute('value')
        self.assertEqual(expected_smiles, smiles_hidden_elem)
        smiles_in_sketcher = sketcher_to_smiles(self)
        self.assertEqual(expected_smiles, smiles_in_sketcher)

    def test_autosave_reaction_table_loaded(self):
        """Tests that the reaction autosaves by reloading the reaction and checking the rxn table data remains"""
        make_new_reaction(self)
        demo_reaction(self)
        # we need to fill in the mass of the limiting reactant
        self.driver.find_element_by_id("js-reactant-rounded-mass1").send_keys(100)
        self.driver.find_element_by_id("js-reactant-equivalent2").send_keys(2)
        time.sleep(1)
        self.exit_and_reload_reaction()
        self.assert_reactants_data()
        self.assert_smiles_data()

    def test_autosave_summary_loaded(self):
        """Tests that the reaction autosaves by editing the summary table and then reloading and checking data remains"""
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        load_status = self.driver.find_element_by_id('js-load-status').get_attribute('value')
        time.sleep(0.5)
        self.assertEqual('loading', load_status)
        time.sleep(6)
        load_status = self.driver.find_element_by_id('js-load-status').get_attribute('value')
        self.assertEqual('loaded', load_status)
        clear_and_send_keys(self, 'js-temperature', '150')
        self.driver.find_element_by_id('cyanide').click()
        time.sleep(1)
        self.exit_and_reload_reaction(reload_id='nitro reduction2-reload')
        # check changed values have autosaved
        cyanide = self.driver.find_element_by_id('cyanide').is_selected()
        self.assertEqual(True, cyanide)
        temperature = self.driver.find_element_by_id('js-temperature').get_attribute('value')
        self.assertEqual('150', temperature)

    def test_autosave_sketcher(self):
        """Tests that the reaction autosaves changes made to the sketcher"""
        make_new_reaction(self)
        self.driver.find_element_by_id('demo-button').click()
        time.sleep(3)
        self.exit_and_reload_reaction()
        time.sleep(6)
        self.assert_smiles_data()

    def test_autosave_after_adding_novel_compound_via_sketcher(self):
        # fill in sketcher and submit with novel compound
        make_new_reaction(self)
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2.CCN>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        self.driver.find_element_by_id('js-new-compound-name').send_keys('TestPhos')
        self.driver.find_element_by_id('js-new-compound-hazards').send_keys('H301')
        self.driver.find_element_by_id('js-new-compound-mw').send_keys('600')
        self.driver.find_element_by_id('js-new-compound-cas').send_keys('123-45-6')
        self.driver.find_element_by_id('js-new-compound-density').send_keys('0.8')
        self.driver.find_element_by_id('js-new-compound-concentration').send_keys('4.0')
        self.driver.find_element_by_id('js-new-compound-submit').click()
        time.sleep(1)
        self.driver.switch_to.alert.accept()
        time.sleep(6)
        # fill reaction table and test that data saves and remains upon reload
        self.driver.find_element_by_id("js-reactant-rounded-mass1").send_keys(100)
        self.driver.find_element_by_id("js-reactant-equivalent2").send_keys(2)
        self.exit_and_reload_reaction()
        self.assert_reactants_data()
        self.assert_smiles_data('O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2.CCN>>C1=CC=CC=C1')

    def test_autosave_reload_with_unfilled_reagent(self):
        """Test autosave works when a unfilled reagent is present"""
        make_new_reaction(self)
        demo_reaction(self)
        self.driver.find_element_by_id("js-reactant-rounded-mass1").send_keys(100)
        self.driver.find_element_by_class_name("js-add-reagent").click()
        time.sleep(1)
        self.exit_and_reload_reaction()
        self.driver.find_element_by_id("js-reactant-equivalent2").send_keys(2)
        self.exit_and_reload_reaction()
        self.assert_reactants_data()
        self.assert_smiles_data()

    def test_autosave_reload_with_unfilled_solvent(self):
        """Test autosave works when a unfilled solvent is present"""
        make_new_reaction(self)
        demo_reaction(self)
        self.driver.find_element_by_id("js-reactant-rounded-mass1").send_keys(100)
        self.driver.find_element_by_class_name("js-add-solvent").click()
        time.sleep(1)
        self.exit_and_reload_reaction()
        self.driver.find_element_by_id("js-reactant-equivalent2").send_keys(2)
        self.exit_and_reload_reaction()
        self.assert_reactants_data()
        self.assert_smiles_data()

    def test_autosave_reload_with_unfilled_add_novel_reagent_to_db(self):
        """Test autosave works when an unfilled novel compound form is present"""
        make_new_reaction(self)
        demo_reaction(self)
        self.driver.find_element_by_id("js-reactant-rounded-mass1").send_keys(100)
        self.driver.find_element_by_id("js-add-new-reagent-by-table").click()
        time.sleep(1)
        self.exit_and_reload_reaction()
        self.driver.find_element_by_id("js-reactant-equivalent2").send_keys(2)
        self.exit_and_reload_reaction()
        self.assert_reactants_data()
        self.assert_smiles_data()

    def test_autosave_after_reloading_sketcher(self):
        """1) make reaction 2) add to sketcher 3) leave 4) reload 5) fill in reaction table 6) reload 7) assert reaction
        table data remains"""
        make_new_reaction(self)
        self.driver.find_element_by_id('demo-button').click()
        time.sleep(3)
        self.exit_and_reload_reaction()
        self.driver.find_element_by_id('action-button-submit').click()
        time.sleep(6)
        self.driver.find_element_by_id("js-reactant-rounded-mass1").send_keys(100)
        self.driver.find_element_by_id("js-reactant-equivalent2").send_keys(2)
        self.exit_and_reload_reaction()
        self.assert_reactants_data()
        self.assert_smiles_data()

    def test_autosave_after_reloading_reaction_table(self):
        make_new_reaction(self)
        demo_reaction(self)
        self.driver.find_element_by_id("js-reactant-rounded-mass1").send_keys(100)
        time.sleep(1)
        self.exit_and_reload_reaction()
        self.driver.find_element_by_id("js-reactant-equivalent2").send_keys(2)
        self.exit_and_reload_reaction()
        self.assert_reactants_data()
        self.assert_smiles_data()

    def test_sketcher_doesnt_autosave_after_reaction_table_formed(self):
        make_new_reaction(self)
        demo_reaction(self)
        rxn_smiles = 'O=C1C2=NC=CN=C2C(=O)C2=C1N=CC=N2>>C1=CC=CC=C1'
        add_reaction_sketcher(self, rxn_smiles, replace=True)
        self.exit_and_reload_reaction()
        self.assert_smiles_data()

    def test_reaction_smiles_update_if_user_submits_second_time(self):
        make_new_reaction(self)
        demo_reaction(self)
        self.driver.find_element_by_id("js-reactant-rounded-mass1").send_keys(100)
        time.sleep(1)
        self.exit_and_reload_reaction()
        rxn_smiles = 'O>>N'
        add_reaction_sketcher(self, rxn_smiles, replace=True)
        time.sleep(1)
        self.driver.find_element_by_id("action-button-submit").click()
        time.sleep(1)
        self.driver.switch_to.alert.accept()
        time.sleep(6)
        self.assert_smiles_data('O>>N')
