from auxiliary_for_tests import *
import flask_testing
from sources import app, db
from database_setup import test_database_create
import unittest
from selenium.webdriver.support.select import Select
from sources import app, db
import time


class ReactTableTestCase(flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        app.config['LOGIN_DISABLED'] = True
        test_database_create(db)
        return app

    def setUp(self):
        pass

    def tearDown(self):
        restore_db()

    def fill_table(self):  # registration function

        return self.client.get(
            '/_process',
            follow_redirects=True
        )

    def filling_table(self, smiles='reactants=C1CCCC1,OC(=O)C1=CC=CC=C1&products=CCNC(=O)C1=CC=CC=C1,N1C=CC=C1'):  # registration function

        return self.client.get(
            '/_process?'+smiles+'&demo="not demo"&'
            'workgroup="Test-Workgroup"&workbook="Test-Workbook"',
            follow_redirects=True
        )

    def send_to_reagents(self, reagent, number):
        return self.client.post(
            '/_reagents',
            data=dict(reagent=reagent, number=number),
            follow_redirects=True
        )

    def test_empty_reaction(self):  # testing successful user registration
        response = self.fill_table()
        self.assertIn(b'{"error":"Missing data!"}\n', response.data)

    def test_reaction_table_product_name(self):  # testing successful user registration
        response = self.filling_table()
        self.assertIn(b'N-Ethylbenzamide', response.data)

    def test_reaction_table_product_mass(self):  # testing successful user registration
        response = self.filling_table()
        self.assertIn(b'value=\\"149.19\\"\\n                       id=\\"js-product-molecular-weight1', response.data)

    def test_reaction_table_reactant_name1(self):  # testing successful user registration
        response = self.filling_table()
        self.assertIn(b'id=\\"js-reactant1\\" style=\\"width: 200px; border-width:0px;\\" value=\\"Cyclopentane\\', response.data)

    def test_reaction_table_reactant_name2(self):  # testing successful user registration
        response = self.filling_table()
        self.assertIn(b'id=\\"js-reactant2\\" style=\\"width: 200px; border-width:0px;\\" value=\\"Benzoic acid\\', response.data)

    def test_reaction_table_reactant_mass1(self):  # testing successful user registration
        response = self.filling_table()
        self.assertIn(b'value=\\"70.13\\"\\n                    id=\\"js-reactant-molecular-weight1', response.data)

    def test_reaction_table_reactant_mass2(self):  # testing successful user registration
        response = self.filling_table()
        self.assertIn(b'value=\\"122.12\\"\\n                    id=\\"js-reactant-molecular-weight2\\', response.data)

    def test_infeasible_reactant(self):
        response=self.filling_table(smiles='reactants=OC(=O)C1=CC=[C](=C)=C=C1,CCN&products=CCNC(=O)C1=CC=CC=C1,N1C=CC=C1')
        self.assertIn(b'Cannot process Reactant 1 structure', response.data)

    def test_infeasible_product(self):
        response=self.filling_table(smiles='reactants=OC(=O)C1=CC=CC=C1,CCN&products=CCNC(=O)C1=CC=CC=C1,CCNC(=O)C1=CC=[C](#C)C=C1')
        self.assertIn(b'Cannot process Product 2 structure', response.data)


class SketcherUnusualStructuresTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    """Tests the reaction table is able to reload when the sketcher contains unusual structures"""

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
        select_workbook(self)
        make_new_reaction(self)

    def test_stereochemistry(self):
        rxn_smiles = "CC[C@@]12CC(C(=O)OC)=C3NC4=CC=CC=C4[C@@]33CCN(CC=C1)[C@@H]23.ClC\C=C\Cl.CC#C>>" \
                     "C[C@H](N)C(O)=O.N[C@@H](Cc1ccccc1)C(O)=O"
        add_reaction_sketcher(self, rxn_smiles)
        self.driver.find_element_by_id("action-button-submit").click()
        time.sleep(8)
        pg_source = self.driver.page_source
        self.assertIn('js-reactant1', pg_source)


if __name__ == '__main__':
    unittest.main()

