from auxiliary_for_tests import *
import flask_testing
from database_setup import test_database_create
from selenium.webdriver.support.select import Select
from sources import app, db
import time


class SketcherReloadTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """ Need to reload the reaction in the sketcher using the data from the reaction database"""

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

    def test_read_marvin_sketch(self):
        """This function tests we can read the contents of the marvin sketcher"""
        select_workgroup(self)
        select_wb = Select(self.driver.find_element_by_id("active-workbook"))
        select_wb.select_by_index(0)
        time.sleep(1)
        make_new_reaction(self)
        self.driver.find_element_by_id("demo-button").click()
        time.sleep(1)
        # function returns the smiles of the structures in the sketcher
        test_smiles = sketcher_to_smiles(self)
        expected_demo_smiles = "OC(=O)C1=CC=CC=C1.CCN>>CCNC(=O)C1=CC=CC=C1"
        self.assertEqual(expected_demo_smiles, test_smiles)

    def test_reload_multiple_reaction_product(self):
        """We test a reaction with multiple reactants and products reloads in the sketcher"""
        select_workgroup(self)
        select_workbook(self, 0)
        time.sleep(1)
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        time.sleep(1)
        # now check reaction has successfully reloaded in the sketcher tool
        test_smiles = sketcher_to_smiles(self)
        expected_smiles = '[O-][N+](=O)C1=CC=CC=C1.O>>NC1=CC=CC=C1.O'
        self.assertEqual(test_smiles, expected_smiles)


class ReactionTableReloadTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """ Need to reload the reaction in the reaction table using the data from the reaction database"""

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
        time.sleep(2)
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        time.sleep(7)

    def test_reactant1_properties(self):
        """Tests reactants properties reload in table"""
        # Reaction table should be loaded now and filled in. Reactant 1 is nitro benzene
        reactant1_name = self.driver.find_element_by_id('js-reactant1').get_attribute('value')
        reactant1_mol_wt = self.driver.find_element_by_id('js-reactant-molecular-weight1').get_attribute('value')
        reactant1_density = self.driver.find_element_by_id('js-reactant-density1').get_attribute('value')
        reactant1_conc = self.driver.find_element_by_id('js-reactant-concentration1').get_attribute('value')
        reactant1_hazards = self.driver.find_element_by_id('js-reactant-hazards1').get_attribute('value')
        self.assertEqual('Nitrobenzene', reactant1_name)
        self.assertEqual('123.11', reactant1_mol_wt)
        self.assertEqual('1.204', reactant1_density)
        self.assertEqual('', reactant1_conc)
        self.assertEqual('H301-H311-H331-H351-H360F-H372-H412', reactant1_hazards)

    def test_reactant2_properties(self):
        # reactant 2 is water
        reactant2_name = self.driver.find_element_by_id('js-reactant2').get_attribute('value')
        reactant2_mol_wt = self.driver.find_element_by_id('js-reactant-molecular-weight2').get_attribute('value')
        reactant2_density = self.driver.find_element_by_id('js-reactant-density2').get_attribute('value')
        reactant2_conc = self.driver.find_element_by_id('js-reactant-concentration2').get_attribute('value')
        reactant2_hazards = self.driver.find_element_by_id('js-reactant-hazards2').get_attribute('value')
        self.assertEqual('Water', reactant2_name)
        self.assertEqual('18.015', reactant2_mol_wt)
        self.assertEqual('0.995', reactant2_density)
        self.assertEqual('', reactant2_conc)
        self.assertEqual('Unknown', reactant2_hazards)

    def test_reactant1_quantities(self):
        """Tests reactant1 quantities reload correctly"""
        reactant1_mass = self.driver.find_element_by_id('js-reactant-rounded-mass1').get_attribute('value')
        reactant1_amount = self.driver.find_element_by_id('js-reactant-rounded-amount1').get_attribute('value')
        reactant1_volume = self.driver.find_element_by_id('js-reactant-rounded-volume1').get_attribute('value')
        self.assertEqual('123', reactant1_mass)
        self.assertEqual('1.00', reactant1_amount)
        self.assertEqual('0.10', reactant1_volume)

    def test_reactant2_quantities(self):
        """Tests reactant2 quantities reload correctly"""
        reactant2_mass = self.driver.find_element_by_id('js-reactant-rounded-mass2').get_attribute('value')
        reactant2_amount = self.driver.find_element_by_id('js-reactant-rounded-amount2').get_attribute('value')
        reactant2_volume = self.driver.find_element_by_id('js-reactant-rounded-volume2').get_attribute('value')
        reactant2_physical_form = self.driver.find_element_by_id('js-reactant-physical-form2').get_attribute('value')
        self.assertEqual('180', reactant2_mass)
        self.assertEqual('10.0', reactant2_amount)
        self.assertEqual('0.02', reactant2_volume)

    def test_reagent1_properties(self):
        # Reaction table should be loaded now and filled in. Reagent1 is iron
        reagent1_name = self.driver.find_element_by_id('js-reagent1').get_attribute('value')
        reagent1_mol_wt = self.driver.find_element_by_id('js-reagent-molecular-weight1').get_attribute('value')
        reagent1_density = self.driver.find_element_by_id('js-reagent-density1').get_attribute('value')
        reagent1_conc = self.driver.find_element_by_id('js-reagent-concentration1').get_attribute('value')
        reagent1_hazards = self.driver.find_element_by_id('js-reagent-hazards1').get_attribute('value')
        self.assertEqual('Methane', reagent1_name)
        self.assertEqual('16.043', reagent1_mol_wt)
        self.assertEqual('', reagent1_density)
        self.assertEqual('', reagent1_conc)
        self.assertEqual('H220', reagent1_hazards)

    def test_reagent2_properties(self):
        # Reaction table should be loaded now and filled in. Reagent2 is Octane
        reagent2_name = self.driver.find_element_by_id('js-reagent2').get_attribute('value')
        reagent2_mol_wt = self.driver.find_element_by_id('js-reagent-molecular-weight2').get_attribute('value')
        reagent2_density = self.driver.find_element_by_id('js-reagent-density2').get_attribute('value')
        reagent2_conc = self.driver.find_element_by_id('js-reagent-concentration2').get_attribute('value')
        reagent2_hazards = self.driver.find_element_by_id('js-reagent-hazards2').get_attribute('value')
        self.assertEqual('Octane', reagent2_name)
        self.assertEqual('114.23', reagent2_mol_wt)
        self.assertEqual('0.703', reagent2_density)
        self.assertEqual('', reagent2_conc)
        self.assertEqual('H225-H302-H304-H312-H315-H319-H336-H400-H410', reagent2_hazards)

    def test_reagent1_quantities(self):
        reagent1_mass = self.driver.find_element_by_id('js-reagent-rounded-mass1').get_attribute('value')
        reagent1_amount = self.driver.find_element_by_id('js-reagent-rounded-amount1').get_attribute('value')
        reagent1_volume = self.driver.find_element_by_id('js-reagent-rounded-volume1').get_attribute('value')
        reagent1_physical_form = self.driver.find_element_by_id('js-reagent-physical-form1').get_attribute('value')
        self.assertEqual('3.21', reagent1_mass)
        self.assertEqual('0.20', reagent1_amount)
        self.assertEqual('', reagent1_volume)
        # self.assertEqual('Dense solid', reagent1_physical_form)

    def test_reagent2_quantities(self):
        reagent2_mass = self.driver.find_element_by_id('js-reagent-rounded-mass2').get_attribute('value')
        reagent2_amount = self.driver.find_element_by_id('js-reagent-rounded-amount2').get_attribute('value')
        reagent2_volume = self.driver.find_element_by_id('js-reagent-rounded-volume2').get_attribute('value')
        reagent2_physical_form = self.driver.find_element_by_id('js-reagent-physical-form2').get_attribute('value')
        self.assertEqual('297', reagent2_mass)
        self.assertEqual('2.60', reagent2_amount)
        self.assertEqual('0.42', reagent2_volume)
        # self.assertEqual('Dense solid', reagent2_physical_form)

    def test_solvent1(self):
        """Tests all the solvent properties are present"""
        # Reaction table should be loaded now and filled in. Solvent 1 is ethanol
        solvent1_name = self.driver.find_element_by_id('js-solvent1').get_attribute('value')
        solvent1_colour = self.driver.find_element_by_id('js-solvent1').value_of_css_property('background-color')
        solvent1_volume = self.driver.find_element_by_id('js-solvent-volume1').get_attribute('value')
        solvent1_conc = self.driver.find_element_by_id('js-solvent-concentration1').get_attribute('value')
        solvent1_hazards = self.driver.find_element_by_id('js-solvent-hazards1').get_attribute('value')
        solvent1_physical_form = self.driver.find_element_by_id('js-solvent-physical-form1').get_attribute('value')
        self.assertEqual('Ethanol', solvent1_name)
        self.assertEqual('2', solvent1_volume)
        self.assertEqual("rgba(0, 255, 0, 1)", solvent1_colour)
        self.assertEqual('H225-H302-H319-H371', solvent1_hazards)
        # self.assertEqual('Highly volatile liquid (b.p. ≤ 70 °C', solvent1_physical_form)

    def test_solvent2(self):
        """Tests all the solvent properties are present"""
        # Reaction table should be loaded now and filled in. Solvent 2 is tetrahydrofuran
        solvent2_name = self.driver.find_element_by_id('js-solvent2').get_attribute('value')
        solvent2_colour = self.driver.find_element_by_id('js-solvent2').value_of_css_property('background-color')
        solvent2_volume = self.driver.find_element_by_id('js-solvent-volume2').get_attribute('value')
        solvent2_conc = self.driver.find_element_by_id('js-solvent-concentration2').get_attribute('value')
        solvent2_hazards = self.driver.find_element_by_id('js-solvent-hazards2').get_attribute('value')
        solvent2_physical_form = self.driver.find_element_by_id('js-solvent-physical-form2').get_attribute('value')
        self.assertEqual('Tetrahydrofuran', solvent2_name)
        self.assertEqual('1', solvent2_volume)
        self.assertEqual(solvent2_colour, "rgba(255, 255, 0, 1)")  # Should be yellow
        self.assertEqual(solvent2_hazards, 'H225-H302-H319-H335-H351')
        # self.assertEqual(solvent2_physical_form, 'Highly volatile liquid (b.p. ≤ 70 °C')

    def test_desired_product(self):
        """Tests the desired product is correct"""
        main_prod = self.driver.find_element_by_id('js-main-product-table-number').get_attribute('value')
        self.assertEqual(main_prod, '1')

    def test_product1(self):
        """Tests all the product properties are present. Product1 is aniline"""
        product_name1 = self.driver.find_element_by_id('js-product1').get_attribute('value')
        product_mol_weight1 = self.driver.find_element_by_id('js-product-molecular-weight1').get_attribute('value')
        product_amount1 = self.driver.find_element_by_id('js-product-rounded-amount1').get_attribute('value')
        product_mass1 = self.driver.find_element_by_id('js-product-rounded-mass1').get_attribute('value')
        product_form1 = self.driver.find_element_by_id('js-product-physical-form1').get_attribute('value')
        product_hazards1 = self.driver.find_element_by_id('js-product-hazard1').get_attribute('value')
        self.assertEqual('Aniline', product_name1)
        self.assertEqual('93.13', product_mol_weight1)
        self.assertEqual('1.00', product_amount1)
        self.assertEqual('93.0', product_mass1)
        self.assertEqual('Gas', product_form1)
        self.assertEqual('H301-H311-H317-H318-H331-H341-H351-H372-H400', product_hazards1)

    def test_product2(self):
        """Tests all the product properties are present. Product2 is water"""
        product_name2 = self.driver.find_element_by_id('js-product2').get_attribute('value')
        product_mol_weight2 = self.driver.find_element_by_id('js-product-molecular-weight2').get_attribute('value')
        product_amount2 = self.driver.find_element_by_id('js-product-rounded-amount2').get_attribute('value')
        product_mass2 = self.driver.find_element_by_id('js-product-rounded-mass2').get_attribute('value')
        product_form2 = self.driver.find_element_by_id('js-product-physical-form2').get_attribute('value')
        product_hazards2 = self.driver.find_element_by_id('js-product-hazard2').get_attribute('value')
        self.assertEqual('Water', product_name2)
        self.assertEqual('18.015', product_mol_weight2)
        self.assertEqual('1.00', product_amount2)
        self.assertEqual('18.0', product_mass2)
        self.assertEqual('Volatile liquid', product_form2)
        self.assertEqual('Unknown', product_hazards2)


class SummaryTableReloadTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """ Need to reload the reaction in the reaction table using the data from the reaction database"""

    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        """We load a test instance of the app, clear and re-populate the database then call the login function"""
        test_database_create(db)
        self.login()

    def tearDown(self):
        # close the browser window
        restore_db()

    def login(self):
        """This function sets the headless browser up, logs the user in and is called from setUp"""
        setup_selenium(self)
        login(self, "BB_Test", "BB_login")
        select_workgroup(self)
        select_workbook(self, 0)
        time.sleep(1)
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        time.sleep(4)

    def value_and_colour(self, cell_id):
        value = self.driver.find_element_by_id(cell_id).get_attribute('value')
        colour = self.driver.find_element_by_id(cell_id).value_of_css_property('background_color')
        return value, colour

    def test_hidden_fields(self):
        """To test the summary table has loaded, the hidden fields are checked"""
        limiting_reactant_mass = self.driver.find_element_by_id('js-limited-reactant-mass').get_attribute('value')
        main_product_mass = self.driver.find_element_by_id('js-main-product-mass').get_attribute('value')
        self.assertEqual("123", limiting_reactant_mass)
        self.assertEqual("93.0", main_product_mass)


class MultipleReloadTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        time.sleep(1)
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        time.sleep(8)

    def test_five_reloads(self):
        """Test reaction preserves data after reloading 5 times"""
        limiting_reactant_mass = self.driver.find_element_by_id('js-limited-reactant-mass').get_attribute('value')
        main_product_mass = self.driver.find_element_by_id('js-main-product-mass').get_attribute('value')
        self.assertEqual("123", limiting_reactant_mass)
        self.assertEqual("93.0", main_product_mass)
        self.check_chem21_metrics()
        # go home and reload reaction
        for i in range(5):
            self.driver.find_element_by_id('TopNavHomeButton').click()
            time.sleep(1)
            select_wg = Select(self.driver.find_element_by_id("WG-select"))
            select_wg.select_by_index(1)
            self.driver.find_element_by_id("go-to-workgroup").click()
            time.sleep(1)
            select_wb = Select(self.driver.find_element_by_id("active-workbook"))
            select_wb.select_by_index(0)
            time.sleep(1)
            self.driver.find_element_by_id('nitro reduction2-reload').click()
            time.sleep(4)
            self.check_chem21_metrics()

    def check_chem21_metrics(self):
        # test chem21 metrics
        temperature = self.driver.find_element_by_id('js-temperature').get_attribute('value')
        elements = self.driver.find_element_by_id('js-elements').get_attribute('value')
        batch_flow = self.driver.find_element_by_id('js-batch-flow').get_attribute('value')
        isolation = self.driver.find_element_by_id('js-isolation').get_attribute('value')
        catalyst = self.driver.find_element_by_id('js-catalyst').get_attribute('value')
        recovery = self.driver.find_element_by_id('js-recovery').get_attribute('value')
        ae = self.driver.find_element_by_id('js-ae-cell').get_attribute('value')
        yield_val = self.driver.find_element_by_id('js-temperature').get_attribute('value')
        self.assertEqual("70", temperature)
        self.assertEqual("+500 years", elements)
        self.assertEqual("Batch", batch_flow)
        self.assertEqual("HPLC", isolation)
        self.assertEqual("Catalyst or enzyme", catalyst)
        self.assertEqual('Not recovered catalyst', recovery)
        self.assertEqual('70', yield_val)
