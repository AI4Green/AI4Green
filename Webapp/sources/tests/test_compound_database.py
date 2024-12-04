from unittest import TestCase, mock, main
import pubchem_download as pcd
import Compound_database_extraction as CDE
from dateutil.relativedelta import relativedelta
from datetime import datetime
from pony.orm import Database, select, db_session, commit, TransactionIntegrityError
import os

"""The first portion of tests here will focus on testing functions inside the Compound_database_extraction module"""


class CompoundDefineEntitiesTestCase(TestCase):
    """Here we are testing the define_database function"""

    def test_add_compound(self):
        """Here we test we can add a compound to the database"""
        db_path = os.path.join('test_data', 'test_db.sqlite')
        self.test_db = Database()
        CDE.define_database(self.test_db)
        self.test_db.bind(provider='sqlite', filename=db_path, create_db=True)
        self.test_db.generate_mapping(create_tables=True)
        with db_session:
            comp1 = self.test_db.Compound(CID=1,
                                          cas='123-45-6',
                                          name='Chemical',
                                          smiles='C=CC',
                                          InChI='InChI=1S/CH4O/c1-2/h2H,1H3',
                                          InChIKey='UHOVQNZJYSORNB-UHFFFAOYSA-N',
                                          molec_formula='C2H5N1',
                                          density=1.07,
                                          boiling_point=110.7,
                                          melting_point=-70.3,
                                          flash_point=300.1,
                                          autoignition_temp=500.8,
                                          molec_weight=112.5,
                                          state='Solid',
                                          form='Crystalline',
                                          hphrase='H123-H456',  # in future this should be an optional StrArray
                                          safety_score=4.1,
                                          health_score=5.2,
                                          enviro_score=1.2,
                                          econom_score=0.2)
            commit()
            a = select(c.name for c in self.test_db.Compound)[:]
            self.assertEqual(a, ['Chemical'])

    def test_wrong_data_type(self):
        """
        Here we provide data that is the wrong type and check an assert error is raised
        This will have to be updated if the database structure is changed.
        """
        db_path = os.path.join('test_data', 'test_db.sqlite')
        self.test_db = Database()
        CDE.define_database(self.test_db)
        self.test_db.bind(provider='sqlite', filename=db_path, create_db=True)
        self.test_db.generate_mapping(create_tables=True)



        """Testing CID only accepts integers."""
        with db_session:
            with self.assertRaises(ValueError):
                comp2 = self.test_db.Compound(CID='not a number',
                                              cas='123-45-6',
                                              name=6)           # This should raise an error, should be a string
                commit()
            """Testing cas code only accepts strings."""
            with self.assertRaises(TypeError):
                comp2 = self.test_db.Compound(CID=1,
                                              cas=123456,
                                              name='chemical')  # This should raise an error, should be a string
                commit()


            """Testing name only accepts strings."""
            with self.assertRaises(TypeError):
                comp2 = self.test_db.Compound(CID=3,
                                              cas='123-45-654',
                                              name=17)  # This should raise an error, should be a string
                commit()

            """Testing smiles only accepts string"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              smiles=90)
                commit()

            """Testing InChI only accepts string"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              InChI=90)
                commit()
            """Testing InChIKey only accepts string"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              InChIKey=['1', 'abc'])
                commit()
            """Testing molecular formula only accepts string"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              molec_formula=56.53)
                commit()
            """Testing density only accepts floats"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              density=[1.234, 135.3])
                commit()

            """Testing boiling point only accepts floats"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              boiling_point='not a float')
                commit()
            """Testing melting point only accepts floats"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              melting_point=['one', 'two'])
                commit()
            """Testing flash point only accepts floats"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              flash_point=['one', 'two'])
                commit()
            """Testing autoignition_temp only accepts floats"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              autoignition_temp=['one', 'two'])
                commit()
            """Testing molec_weight only accepts floats"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              molec_weight=['one', 'two'])
                commit()
            """Testing hazards only accepts strings"""
            with self.assertRaises(TypeError):
                comp1 = self.test_db.Compound(CID=2,
                                              cas='123-45-689',
                                              name='Chemical',
                                              hphrase=8)
                commit()

    def test_unique_compound(self):
        """Here we test that the compound database won't accept non-unique entries for CID and CAS"""
        db_path = os.path.join('test_data', 'test_db.sqlite')
        self.test_db = Database()
        CDE.define_database(self.test_db)
        self.test_db.bind(provider='sqlite', filename=db_path, create_db=True)
        self.test_db.generate_mapping(create_tables=True)
        with db_session:
            comp1 = self.test_db.Compound(CID='1',
                                          cas='123-45-6',
                                          name='Chemical',
                                          smiles='C=CC',
                                          InChI='InChI=1S/CH4O/c1-2/h2H,1H3',
                                          InChIKey='UHOVQNZJYSORNB-UHFFFAOYSA-N',
                                          molec_formula='C2H5N1',
                                          density=1.07,
                                          boiling_point=110.7,
                                          melting_point=-70.3,
                                          flash_point=300.1,
                                          autoignition_temp=500.8,
                                          molec_weight=112.5,
                                          state='Solid',
                                          form='Crystalline',
                                          hphrase='H123-H456',  # in future this should be an optional StrArray
                                          safety_score=4.1,
                                          health_score=5.2,
                                          enviro_score=1.2,
                                          econom_score=0.2)
            commit()

        """Now we add another compound with the same CAS and we assert an error from pony.orm is raised"""

        with self.assertRaises(TransactionIntegrityError):
            with db_session:
                comp2 = self.test_db.Compound(CID='2',
                                              cas='123-45-6',
                                              name='Chemical',
                                              smiles='C=CC',
                                              InChI='InChI=1S/CH4O/c1-2/h2H,1H3',
                                              InChIKey='UHOVQNZJYSORNB-UHFFFAOYSA-N',
                                              molec_formula='C2H5N1',
                                              density=1.07,
                                              boiling_point=110.7,
                                              melting_point=-70.3,
                                              flash_point=300.1,
                                              autoignition_temp=500.8,
                                              molec_weight=112.5,
                                              state='Solid',

                                              form='Crystalline',
                                              hphrase='H123-H456',  # in future this should be an optional StrArray
                                              safety_score=4.1,
                                              health_score=5.2,
                                              enviro_score=1.2,
                                              econom_score=0.2)
                commit()

        with self.assertRaises(TransactionIntegrityError):
            with db_session:
                comp2 = self.test_db.Compound(CID='1',
                                              cas='123-45-90183',
                                              name='Chemical',
                                              smiles='C=CC',
                                              InChI='InChI=1S/CH4O/c1-2/h2H,1H3',
                                              InChIKey='UHOVQNZJYSORNB-UHFFFAOYSA-N',
                                              molec_formula='C2H5N1',
                                              density=1.07,
                                              boiling_point=110.7,
                                              melting_point=-70.3,
                                              flash_point=300.1,
                                              autoignition_temp=500.8,
                                              molec_weight=112.5,
                                              state='Solid',

                                              form='Crystalline',
                                              hphrase='H123-H456',  # in future this should be an optional StrArray
                                              safety_score=4.1,
                                              health_score=5.2,
                                              enviro_score=1.2,
                                              econom_score=0.2)
                commit()


    def tearDown(self):
        self.test_db.disconnect()
        outfile1 = os.path.join('test_data/test_db.sqlite')
        os.remove(outfile1)


class CompoundDatabaseSmilesTestCase(TestCase):
    """Here we are testing the function that turns inchi into SMILES"""

    def test_conversion_methanol(self):
        test_inchi = 'InChI=1S/CH4O/c1-2/h2H,1H3'
        test_smiles_str = CDE.get_smiles(test_inchi)
        known_smiles_str = 'CO'
        self.assertEqual(test_smiles_str, known_smiles_str)

    # Seems to rely on chemspider - limited API usage so is off by default
    # def test_conversion_complex_aromatic(self):
    #     test_inchi = 'InChI=1S/C25H24N6O2/c1-2-21(32)30-14-6-7-18(15-30)31-25-22(24(26)27-16-28-25)23(' \
    #                  '29-31)17-10-12-20(13-11-17)33-19-8-4-3-5-9-19/h2-5,8-13,16,18H,1,6-7,14-15H2,(H2,26,27,' \
    #                  '28)/t18-/m1/s1 '
    #     test_smiles_str = CDE.get_smiles(test_inchi)
    #     known_smiles_str = 'C=CC(=O)N1CCC[C@H](C1)n1c2c(c(c3ccc(cc3)Oc3ccccc3)n1)c(N)ncn2'
    #     self.assertEqual(test_smiles_str, known_smiles_str)

    def test_conversion_caffeine(self):
        test_inchi = 'InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3'
        test_smiles_str = CDE.get_smiles(test_inchi)
        known_smiles_str = 'Cn1c(=O)c2c(ncn2C)n(C)c1=O'
        self.assertEqual(test_smiles_str, known_smiles_str)



class CompoundDatabaseExtractionTestCase(TestCase):
    """Here we are testing the compound database extraction function"""

    def test_extraction(self):
        """
        In this function we want to test the function is capable of pulling out all of the attributes from the pubchem
        style xml and saving them to the database.
        """

        self.test_db = Database()
        CDE.define_database(self.test_db)
        self.test_db.bind(provider='sqlite', filename='test_data/test_db.sqlite', create_db=True)
        self.test_db.generate_mapping(create_tables=True)
        CDE.extract_from_pubchem('test_data/Old_haz_data.xml', self.test_db, 10)
        with db_session:
            a = select((c.CID, c.id, c.hphrase, c.cas, c.smiles, c.InChI, c.name,
                        c.molec_formula, c.molec_weight) for c in self.test_db.Compound)[:]
            for attributes in a:
                self.assertIsNotNone(attributes)

    def tearDown(self):
        self.test_db.disconnect()
        outfile1 = os.path.join('test_data/test_db.sqlite')
        os.remove(outfile1)


class DatabaseUpdateStage1TestCase(TestCase):
    """Here we test the code for stage1 of the update. This function is called then the database status
    is set to 'Update'. Normally we use the pubchem url to retrieve a list of update and return all of
    the updates which have been released since we last updated our database"""

    @mock.patch('pubchem_download.find_files')
    @mock.patch('pubchem_download.get_date')
    def test_stage1(self, mock_date, mock_ff):
        """Inside this function we test the function responds as expected when changing the input"""
        """
        We want to return all of these updates because they are all more recent than the date of our previous update
        """
        print("Testing database update stage 1")
        mock_date.return_value = datetime(2021, 1, 1)
        mock_ff.return_value = ['PubChem-updates', '2021-02-01/', '2021-03-01/', '2021-04-01/']
        result = pcd.find_database_updates()
        self.assertEqual(result, ['2021-02-01/', '2021-03-01/', '2021-04-01/'], "Wrong list obtained")

        """We want to test that any older values are not returned, so here we only want the most recent date"""
        mock_ff.return_value = ['PubChem-updates', '2020-02-01/', '2020-03-01/', '2021-04-01/']
        result = pcd.find_database_updates()
        self.assertEqual(result, ['2021-04-01/'], "Wrong list obtained")

        """We want to test that an error is raised if an object which can't be converted to datetime is in the list"""
        mock_ff.return_value = ['PubChem-updates', 'Not a date', 'Also not a date']
        result = pcd.find_database_updates()
        self.assertEqual(result, [])


@mock.patch('pubchem_download.md5_integrity_check')
@mock.patch('pubchem_download.download_file')
@mock.patch('pubchem_download.unzip_file')
@mock.patch('pubchem_download.clear_update_table')
class DatabaseUpdateStage2TestCase(TestCase):
    """Here we are testing stage 2
    This function will download the update database and checksum, check the integrity and unzip the database xml,
    extract the xml into an sqlite database and then compare this to the old database and update any differences to the
    old database"""

    def tearDown(self):
        self.test_db3.disconnect()
        self.test_db2.disconnect()
        self.test_db1.disconnect()
        outfile1 = os.path.join('test_data/test_db1.sqlite')
        outfile2 = os.path.join('test_data/test_db2.sqlite')
        outfile3 = os.path.join('test_data/test_db3.sqlite')
        os.remove(outfile1)
        os.remove(outfile2)
        os.remove(outfile3)

    def setUp(self):
        """Here we set up the three test databases"""

        def setup_test_db(in_text_file, out_db_file):
            db = Database()
            CDE.define_database(db)
            db.bind(provider='sqlite', filename=out_db_file, create_db=True)
            db.generate_mapping(create_tables=True)
            CDE.extract_from_pubchem(in_text_file, db, 10)
            print(out_db_file, "set up")
            return db

        """This data is a single compound"""
        old_haz_data = os.path.join('test_data/Old_haz_data.xml')
        outfile1 = os.path.join('test_data/test_db1.sqlite')
        self.test_db1 = setup_test_db(old_haz_data, outfile1)

        """This data is a single compound with updated hazard codes"""
        new_haz_data = os.path.join('test_data/New_haz_data.xml')
        outfile2 = os.path.join('test_data/test_db2.sqlite')
        self.test_db2 = setup_test_db(new_haz_data, outfile2)

        """This data is the single compound from test_db1 plus 1 extra compound"""
        new_cmpd_data = os.path.join('test_data/New_cmpd_data.xml')
        outfile3 = os.path.join('test_data/test_db3.sqlite')
        self.test_db3 = setup_test_db(new_cmpd_data, outfile3)

    @mock.patch('Compound_database_extraction.extract_from_pubchem')
    @mock.patch('pubchem_download.get_date')
    @mock.patch('pubchem_download.read_checksum_file')
    def test_stage2_func_calls(self, mock_checksum_verify, mock_date, mock_extract, mock_clear, mock_unzip, mock_dl, mock_md5):
        """In this function we test all the expected function calls are made after the test_stage2 function is called"""

        self.assertFalse(mock_dl.called, "Must downloaded the database update files")
        self.assertFalse(mock_unzip.called, "Must unzip the downloaded files")
        self.assertFalse(mock_extract.called, "Must extract the data into the update database")
        self.assertFalse(mock_clear.called), "Must clear the table of the update database"

        """We need to setup the rest of the variables"""
        today = datetime.today()

        past_date = today.strptime('2021-01-01', '%Y-%m-%d')
        mock_md5.return_value = 'Pass'

        pcd.apply_database_updates(self.test_db1, self.test_db2)

    @mock.patch('Compound_database_extraction.extract_from_pubchem')
    def test_stage2_update_existing(self, mock_extract, mock_clear, mock_unzip, mock_dl, mock_md5):
        """In this function we want to test updating the hazard codes for an existing entry in the database
        self.test_db1: is the old database
        self.test_db2: is the update database
        """
        print("Testing update stage 2 for updating existing compounds")

        with db_session:
            a = select(c.hphrase for c in self.test_db1.Compound if c.id == 1)[:]
            b = select(c.hphrase for c in self.test_db2.Compound if c.id == 1)[:]
            self.assertNotEqual(['H999-H998-H997-H996-H995-H994-H999'], a, "test_db1 must not have these codes yet")
            self.assertNotEqual(a, b, "The two hazards must be different at the start of testing")

        """We need to setup the rest of the variables"""
        today = datetime.today()
        today = today.strftime('%Y-%m-%d')
        mock_md5.return_value = 'Pass'

        pcd.apply_database_updates(self.test_db1, self.test_db2)

        """We check that these functions are all called as part of database_update_stage2 
        and that the hazard codes of test_db1 are now == to those of test_db2"""

        with db_session:
            a = select(c.hphrase for c in self.test_db1.Compound if c.id == 1)[:]
            self.assertEqual(['H999-H998-H997-H996-H995-H994-H999'], a, "test db1 should be updated to these hazards")
            self.assertEqual(a, b, "The two hazards must be the same after database_update_stage2")

    @mock.patch('Compound_database_extraction.extract_from_pubchem')
    def test_stage2_add_new(self, mock_extract, mock_clear, mock_unzip, mock_dl, mock_md5):
        """In this function we want to add a new compound to the old database
        self.test_db1: is the old database
        self.test_db2 is the update database
        """
        print("Testing update stage 2 for adding new compounds")
        """        
        First we test that the setup has worked as intended and that the compounds in db1 and db3 are different
        """
        with db_session:
            a = select(c.CID for c in self.test_db1.Compound)[:]
            b = select(c.CID for c in self.test_db3.Compound)[:]
            self.assertNotEqual(a, b, "The two databases must have different compounds in at this stage")
            self.assertEqual(len(a), 1, "Test db1 must only have a single entry")
            self.assertEqual(len(b), 2, "Test db3 must have 2 entries")

        """We need to assign the rest of the variables and then run the function"""
        today = datetime.today()
        today = today.strftime('%Y-%m-%d')
        mock_md5.return_value = 'Pass'

        pcd.apply_database_updates(self.test_db1, self.test_db3)

        """Here we test that the new entry has been added to test_db1"""
        with db_session:
            a = select(c.CID for c in self.test_db1.Compound)[:]
            self.assertEqual(a, b, "The two databases must now contain the same compounds")
            self.assertEqual(len(a), 2, "Test db1 must have 2 entries now")

    @mock.patch('Compound_database_extraction.extract_from_pubchem')
    def test_stage2_del_old(self, mock_extract, mock_clear, mock_unzip, mock_dl, mock_md5):
        """In this function we want to test deleting an old compound
        self.test_db3: is the old database
        self.test_db1: is the update database
        """
        print("Testing update stage 2 for deleting compounds")
        """        
        First we test that the setup has worked as intended and that the compounds in db1 and db3 are different
        """
        with db_session:
            a = select(c.CID for c in self.test_db3.Compound)[:]
            b = select(c.CID for c in self.test_db1.Compound)[:]
            self.assertNotEqual(a, b, "The two databases must have different compounds in at this stage")
            self.assertEqual(len(a), 2, "Test db3 must have 2 entries")
            self.assertEqual(len(b), 1, "Test db1 must only have 1 entry")

        """We need to assign the rest of the variables and then run the function"""
        today = datetime.today()
        today = today.strftime('%Y-%m-%d')
        mock_md5.return_value = 'Pass'

        pcd.apply_database_updates(self.test_db3, self.test_db1)

        """Here we test the old entry has been deleted from test_db3"""
        with db_session:
            a = select(c.CID for c in self.test_db3.Compound)[:]
            self.assertEqual(a, b, "The two databases must now contain the same compounds")
            self.assertEqual(len(a), 1, "Test db3 must only have 1 entries now")




if __name__ == '__main__':
    main()
