from unittest import TestCase, mock, main
from pony.orm import commit, db_session, select, set_sql_debug, Database, TransactionIntegrityError
import Reactions_db_define as RXNBD
import logging
import sys
import re
import pathlib as pl
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta as timediff

regex = '^[a-z0-9]+[\._]?[a-z0-9]+[@]\w+[.]\w{2,3}$'




class lookforDatabase(TestCase):
    """
    This tests that the database is built and present in the first place
    """
    def test_build_database(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        path = pl.Path(now + "sqlite")
        self.assertFalse(os.path.isfile(path))

    def tearDown(self):
        pass
        self.db_test.drop_all_tables(with_all_data=True)



class addNovelCompoundTest(TestCase):
    """
    This class is to test the adding of a Novel compound to the Reaction database

    In future this test must be expanded so that;

    a) a User can only add a novel compound to their workbooks that they have access to
    b) a User can only access a novel compound in their workbook libraries
    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(True)
        #log = logging.getLogger("SomeTest.testSomething")
        #log.debug("test")
        #log.debug("testtest")
        add_users(self.db_test)

    @db_session
    def test_add_compound_1(self):
        '''
        This test looks at adding one novel compound with the required attributes
        name
        group

        The other attributes are not required for a novel compound
        '''
        wg1 = select(c.workbooks for c in self.db_test.WorkGroup if c.name == "Smith_Lab")[:]
        wg1 = wg1[0]
        wg2 = select(c.workbooks for c in self.db_test.WorkGroup if c.name == "Willis_Lab")[:]
        wg2 = wg2[0]

        nc1 = self.db_test.NovelCompound(name="PixieDust", workbook=wg1)
        nc2 = self.db_test.NovelCompound(name="Unobtainium", workbook=wg2)


        commit()

        result = select(c for c in self.db_test.NovelCompound if c.name == "PixieDust")[:]
        resultname = result[0].name
        resultworkbook = result[0].workbook
        self.assertEqual(resultname, "PixieDust")
        self.assertEqual(resultworkbook, wg1)
        result = select(c for c in self.db_test.NovelCompound if c.name == "Unobtainium")[:]
        resultname = result[0].name
        resultworkbook = result[0].workbook
        self.assertEqual(resultname, "Unobtainium")
        self.assertEqual(resultworkbook, wg2)


    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)

class addNonUniqueNovelCompoundTest(TestCase):
    """
    This class is to test the adding of a Novel compound to the Reaction database

    In future this test must be expanded so that;

    a) a User can only add a novel compound to their workbooks that they have access to
    b) a User can only access a novel compound in their workbook libraries
    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(True)
        #log = logging.getLogger("SomeTest.testSomething")
        #log.debug("test")
        #log.debug("testtest")
        add_users(self.db_test)


    def test_add_compound_1(self):
        with db_session:
            wg1 = select(c.workbooks for c in self.db_test.WorkGroup if c.name == "Smith_Lab")[:]
            wg1 = wg1[0]
            nc1 = self.db_test.NovelCompound(name="PixieDust", workbook=wg1)

        with self.assertRaises(TransactionIntegrityError):
            with db_session:

                '''
                This test looks at adding one non unique novel compound with the required attributes
                name
                group
    
                The other attributes are not required for a novel compound
    
                It should throw an error
                '''
                wg1 = select(c.workbooks for c in self.db_test.WorkGroup if c.name == "Smith_Lab")[:]
                wg1 = wg1[0]
                nc3 = self.db_test.NovelCompound(name="PixieDust", workbook=wg1)

                commit()

    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)

class addNovelCompound_Full_Test(TestCase):
    """
    This class is to test the adding of a Novel compound to the Reaction database

    This tests that we can fill all the data attributes
    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(True)
        add_users(self.db_test)


    def test_add_compound_1(self):
        with self.assertRaises(TypeError):
            with db_session:
                wg1 = select(c.workbooks for c in self.db_test.WorkGroup if c.name == "Smith_Lab")[:]
                wg1 = wg1[0]
                nc1 = self.db_test.NovelCompound(name="PixieDust", workbook=wg1, CAS="1223",
                                                 molec_formula=1.3, density="1.3", boiling_point="100.1",
                                                 melting_point="20.1", flash_point="Test", autoignition_temp="70.2",
                                                 molec_weight="25.0", state=0.5, form=22.3, Hazards=[23.4, 56],
                                                 safety_score="23.4", health_score="5.5", enviro_score="6.7",
                                                 econom_score="4.3")
                commit()



    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)


class addNewUser(TestCase):
    """
    This class is to test the adding of a user to the user Compound database

    In future this test must be expanded so that;

    a) a User can be added only by a PI
    b) a User Must be associated with a Workgroup
    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(True)
        add_users(self.db_test)

    @db_session
    def test_add_user(self):
        '''
        This test looks at adding one user

        '''
        wg1 = select(c for c in self.db_test.WorkGroup if c.name == "Smith_Lab")[:]
        wg1 = wg1[0]
        wg2 = select(c for c in self.db_test.WorkGroup if c.name == "Willis_Lab")[:]
        wg2 = wg2[0]

        p1 = self.db_test.Person(firstname='Peter', middlename='Paul', lastname='Parker', email='ppparker@gmail.com')


        commit()

        user1 = select(c for c in self.db_test.Person if c.lastname == "Parker")[:]
        user1 = user1[0]
        self.assertEqual(user1.firstname, "Peter")
        self.assertEqual(user1.middlename, "Paul")
        self.assertEqual(user1.lastname, "Parker")
        user1email = check(user1.email)
        self.assertEqual(True, user1email)

    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)

class addNonUniqueNewUser(TestCase):
    """
    This class is to test the adding of a user to the user Compound database

    In future this test must be expanded so that;

    a) a User can be added only by a PI
    b) a User Must be associated with a Workgroup
    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(True)
        add_users(self.db_test)

    def test_add_user(self):
        with db_session:
            wg1 = select(c for c in self.db_test.WorkGroup if c.name == "Smith_Lab")[:]
            wg1 = wg1[0]
            p1 = self.db_test.Person(firstname='Peter', middlename='Paul', lastname='Parker', email='ppparker@gmail.com')
            commit()
        '''
        This test looks at adding one user

        '''
        with self.assertRaises(TransactionIntegrityError):
            with db_session:
                p2 = self.db_test.Person(firstname='Phil', middlename='Paul', lastname='Parker', email='ppparker@gmail.com')
                commit()

    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)




class addNewWorkgroup(TestCase):
    """
    This class is to test the adding of a Workgroup to the Reactions Database

    In future this test must be expanded so that;

    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(False)
        add_users(self.db_test)

    @db_session
    def test_new_workgroup(self):
        '''
        This test looks at adding one work group
        '''

        useremail = "dsbarlow@gmail.com"
        user1email = check(useremail)

        p1 = select(c for c in self.db_test.Person if c.email == useremail)[:]
        i1 = select(c for c in self.db_test.Institution if c.name == 'Uni Test')[:][0]
        wg1 = self.db_test.WorkGroup(name='Barlow_Lab', principal_investigator=p1, institution=i1)

        commit()


    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)

class addNonUniqueNewWorkgroup(TestCase):
    """
    This class is to test the adding of a non unique Workgroup to the Reactions Database

    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(False)
        add_users(self.db_test)

    def test_new_workgroup(self):
        '''
        This test looks at adding one work group
        '''
        useremail = "brwillis@gmail.com"
        user1email = check(useremail)

        useremail2 = "jjsmith@gmail.com"
        user1email = check(useremail2)

        with db_session:
            p2 = select(c for c in self.db_test.Person if c.email == useremail2).first()

            p22 = select(d for c in self.db_test.Person
                         for d in c.workgroup_principal_investigator
                         if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

            wb4 = self.db_test.WorkBook(name='Project3', group=p22, users=p2)
            commit()
        with self.assertRaises(TransactionIntegrityError):
            with db_session:
                p2 = select(c for c in self.db_test.Person if c.email == useremail2).first()

                p22 = select(d for c in self.db_test.Person
                             for d in c.workgroup_principal_investigator
                             if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

                wb6 = self.db_test.WorkBook(name='Project3', group=p22, users=p2)
                commit()


    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)


class addNewWorkbook(TestCase):
    """
    This class is to test the adding of a Workbook to the Reactions Database

    In future this test must be expanded so that;

    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(False)
        add_users(self.db_test)

    @db_session
    def test_new_workbook(self):
        '''
        This test looks at adding one workbook
        '''

        useremail = "brwillis@gmail.com"
        user1email = check(useremail)

        useremail2 = "jjsmith@gmail.com"
        user1email = check(useremail2)

        p2 = select(c for c in self.db_test.Person if c.email == useremail2).first()

        p22 = select(d for c in self.db_test.Person
                     for d in c.workgroup_principal_investigator
                     if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

        # print(p2)
        # print(p22)

        wb4 = self.db_test.WorkBook(name='Project3', group=p22, users=p2)


        commit()


    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)

class addNewReaction(TestCase):
    """
    This class is to test the adding of a Reaction to the Reactions Database

    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(False)
        add_users(self.db_test)

    @db_session
    def test_new_reaction(self):

        useremail = "brwillis@gmail.com"
        user1email = check(useremail)

        useremail2 = "jjsmith@gmail.com"
        user1email = check(useremail2)

        p2 = select(c for c in self.db_test.Person if c.email == useremail2).first()

        wg1 = select(d for c in self.db_test.Person
                     for d in c.workgroup_principal_investigator
                     if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

        wb = select(e for c in self.db_test.Person
                     for d in c.workgroup_principal_investigator
                     for e in d.workbooks
                     if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

        # print(p2)
        # print(p22)

        wb1 = self.db_test.WorkBook(name='Project3', group=wg1, users=p2)

        wb2 = select(c for c in self.db_test.WorkBook if c.name == 'Project3' and c.group == wg1).first()

        today = datetime.today()
        #print(today)
        futureday = today + timediff(months=1)
        #print(futureday)

        rxn1 = self.db_test.Reaction(name="Test1", time_of_creation=today,
                                     date_reaction=futureday, workbooks=wb1,
                                     reactants=[1,2], products=[5,6], creator_email=p2.email)

        rxn11 = select(c for c in self.db_test.Reaction if c.name == "Test1").first()

        commit()


    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)

class addNonUniqueReaction(TestCase):
    """
        This test looks at adding a non-unique reaction

    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(False)
        add_users(self.db_test)

    def test_new_reaction(self):

        useremail = "brwillis@gmail.com"
        user1email = check(useremail)

        useremail2 = "jjsmith@gmail.com"
        user1email = check(useremail2)

        today = datetime.today()
        futureday = today + timediff(months=1)


        with db_session:
            p2 = select(c for c in self.db_test.Person if c.email == useremail2).first()

            p22 = select(d for c in self.db_test.Person
                         for d in c.workgroup_principal_investigator
                         if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

            p222 = select(e for c in self.db_test.Person
                         for d in c.workgroup_principal_investigator
                         for e in d.workbooks
                         if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

            wb1 = self.db_test.WorkBook(name='Project3', group=p22, users=p2)

            wb2 = select(c for c in self.db_test.WorkBook if c.name == 'Project3' and c.group == p22).first()

            rxn1 = self.db_test.Reaction(name="Test1", time_of_creation=today,
                                        date_reaction=futureday, workbooks=wb1,
                                         reactants=[1, 2], products=[5, 6], creator_email=p2.email)
            commit()


        with self.assertRaises(TransactionIntegrityError):
            with db_session:
                p2 = select(c for c in self.db_test.Person if c.email == useremail2).first()

                p22 = select(d for c in self.db_test.Person
                             for d in c.workgroup_principal_investigator
                             if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

                p222 = select(e for c in self.db_test.Person
                              for d in c.workgroup_principal_investigator
                              for e in d.workbooks
                              if (c.email == useremail2) and (d.name == "Smith_Lab")).first()


                wb1 = select(c for c in self.db_test.WorkBook
                             if (c.group == p22) and (c.name == "Project3")).first()

                rxn2 = self.db_test.Reaction(name="Test1", time_of_creation=today,
                                             date_reaction=futureday, workbooks=wb1,
                                             reactants=[4, 3], products=[5, 6], creator_email=p2.email)

                commit()



    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)

class addLinkedReaction(TestCase):
    """
    This class is to test the adding of a Reaction linked to another Reaction to the Reactions Database

    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(True)
        add_users(self.db_test)

    @db_session
    def test_new_reaction(self):

        useremail = "brwillis@gmail.com"
        user1email = check(useremail)

        useremail2 = "jjsmith@gmail.com"
        user1email = check(useremail2)

        p2 = select(c for c in self.db_test.Person if c.email == useremail2).first()

        p22 = select(d for c in self.db_test.Person
                     for d in c.workgroup_principal_investigator
                     if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

        p222 = select(e for c in self.db_test.Person
                     for d in c.workgroup_principal_investigator
                     for e in d.workbooks
                     if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

        # print(p2)
        # print(p22)

        wb1 = self.db_test.WorkBook(name='Project3', group=p22, users=p2)

        wb2 = select(c for c in self.db_test.WorkBook if c.name == 'Project3' and c.group == p22).first()

        today = datetime.today()
        #print(today)
        futureday = today + timediff(months=1)
        #print(futureday)

        rxn1 = self.db_test.Reaction(name="Test1", time_of_creation=today,
                                     date_reaction=futureday, workbooks=wb1,
                                     reactants=[1,2], products=[5,6], creator_email=p2.email)

        rxn11 = select(c for c in self.db_test.Reaction if c.name == "Test1").first()

        commit()

        today = futureday + timediff(days=5)
        futureday = today + timediff(months=2)

        rxn2 = self.db_test.Reaction(name="Test2", time_of_creation=today, date_reaction=futureday,
                                     creator_email=p2.email, workbooks=wb1, reactants=[6,7], products=[8,9],
                                     precursor_reaction=rxn11)

        rxn1.successor_reaction = rxn2

        commit()



    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)

class addWrongLinkedReaction(TestCase):
    """
    This class is to test the adding of a Reaction to the Reactions Database using the wrong type for linking

    """

    def setUp(self):
        now = 'reactions_database_test'
        self.db_test = RXNBD.open_database(now)
        set_sql_debug(True)
        add_users(self.db_test)

    def test_new_reaction(self):

        useremail = "brwillis@gmail.com"
        user1email = check(useremail)

        useremail2 = "jjsmith@gmail.com"
        user1email = check(useremail2)


        with db_session:
            p2 = select(c for c in self.db_test.Person if c.email == useremail2).first()

            p22 = select(d for c in self.db_test.Person
                         for d in c.workgroup_principal_investigator
                         if (c.email == useremail2) and (d.name == "Smith_Lab")).first()

            p222 = select(e for c in self.db_test.Person
                         for d in c.workgroup_principal_investigator
                         for e in d.workbooks
                         if (c.email == useremail2) and (d.name == "Smith_Lab")).first()


            wb1 = self.db_test.WorkBook(name='Project3', group=p22, users=p2)

            wb2 = select(c for c in self.db_test.WorkBook if c.name == 'Project3' and c.group == p22).first()

            today = datetime.today()
            futureday = today + timediff(months=1)

            rxn1 = self.db_test.Reaction(name="Test1", time_of_creation=today, creator_email=p2.email,
                                         date_reaction=futureday, workbooks=wb1,
                                         reactants=[1,2], products=[5,6])


            commit()

        with self.assertRaises(TypeError):
            with db_session:
                today = futureday + timediff(days=5)
                futureday = today + timediff(months=2)
                rxn11 = select(c for c in self.db_test.Reaction if c.name == "Test1").first()

                rxn2 = self.db_test.Reaction(name="Test2", time_of_creation=today, date_reaction=futureday,
                                             workbooks=wb1, reactants=[6,7], products=[8,9], creator_email=p2.email,
                                             precursor_reaction="a")


                commit()


    def tearDown(self):
        self.db_test.drop_all_tables(with_all_data=True)



@db_session
def add_users(db):
    """Populates the database with users, workgroups, workbooks and forms the connections between these"""

    p1 = db.Person(firstname='John', middlename='Jay', lastname='Smith', email='jjsmith@gmail.com')
    p2 = db.Person(firstname='Mary', middlename='May', lastname='Thompson', email='mmthompson@gmail.com')
    p3 = db.Person(firstname='Bill', middlename='Roger', lastname='Willis', email='brwillis@gmail.com')
    p4 = db.Person(firstname='Donna', middlename='Sue', lastname='Barlow', email='dsbarlow@gmail.com')

    i1 = db.Institution(name='Uni Test')

    wg1 = db.WorkGroup(name='Smith_Lab', principal_investigator=p1, standard_member=[p1,p2,p3], institution=i1)
    wg2 = db.WorkGroup(name='Willis_Lab', principal_investigator=p3, standard_member=[p3,p4], institution=i1)


    wb1 = db.WorkBook(name='Project1', group=wg1, users=[p1,p2,p3])
    wb2 = db.WorkBook(name='Project2', group=wg2, users=[p3,p4])
    wb3 = db.WorkBook(name='Project3', group=wg2, users=p3)

    commit()





def check(email):
    if (re.search(regex, email)):
        #print("Valid Email")
        return True
    else:
        #print("Invalid Email")
        return False


if __name__ == '__main__':
    print('Testing')
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("SomeTest.testSomething").setLevel(logging.DEBUG)
    main()
