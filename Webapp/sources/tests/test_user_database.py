from auxiliary_for_tests import *
from pony.orm import select, commit, db_session
from pony.orm.core import TransactionIntegrityError
from sources import app, db
import re
from werkzeug.security import generate_password_hash, check_password_hash
from database_setup import test_database_create
import unittest

regex = '^[a-z0-9]+[\._]?[a-z0-9]+[@]\w+[.]\w{2,3}$'


class UserTestCase(unittest.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        test_database_create(db)

    def tearDown(self):
        restore_db()

    @db_session
    def test_check_user(self):  # checks if the user is in the database
        user = db.User.get(username='BB_Test')
        self.assertEqual(user.email, 'BB@test.com')

    @db_session
    def test_add_unique_user(self):  # checks if a unique user can be added in the database
        p1 = db.Person()
        user = db.User(fullname='Thomas Thomson', username='Tom', email='tom@mail.com',
                       person=p1, password_hash=generate_password_hash('pass'))
        a = select(u.email for u in db.User)[:]
        self.assertEqual(['BB@test.com', 'PI@test.com', "SM1@test.com", 'SM@test.com', 'SR@test.com', 'admin@test.com',
                          'tom@mail.com'], a)

    @db_session
    def test_nonunique_username(self):  # checks if a user with a nonunique username can be added in the database
        with self.assertRaises(TransactionIntegrityError):
            p1 = db.Person()
            db.User(fullname='Bob Brown', username='BB_Test', email='joe@mail.com',
                    person=p1, password_hash=generate_password_hash('hash'))
            commit()

    @db_session
    def test_nonunique_email(self):  # checks if a user with a nonunique email can be added in the database
        with self.assertRaises(TransactionIntegrityError):
            p1 = db.Person()
            db.User(fullname='Joseph Jones', username='Joe', email='BB@test.com',
                    person=p1, password_hash=generate_password_hash('hash'))
            commit()

    @db_session
    def test_good_email(self):  # checks if the email is in the correct format
        p1 = db.Person()
        user = db.User(fullname='Thomas Thomson', username='Tom', email='tom@mail.com',
                       person=p1, password_hash=generate_password_hash('pass'))
        user = db.User.get(username='Tom')
        self.assertEqual(check(user.email), True)

    @db_session
    def test_bad_email(self):  # checks if an email is in an incorrect format
        p1 = db.Person()
        user = db.User(fullname='Thomas Thomson', username='Tom', email='tom@mail',
                       person=p1, password_hash=generate_password_hash('pass'))
        self.assertEqual(check(user.email), False)

    @db_session
    def test_check_password(self):  # checks if the password is hashed
        user = select(u for u in db.User if 'BB_Test' == u.username).first()
        a = check_password_hash(user.password_hash, 'BB_login')
        self.assertEqual(a, True)


def check(email):
    if re.search(regex, email):
        return True
    else:
        return False


if __name__ == '__main__':
    unittest.main()
