from auxiliary_for_tests import *
from flask_testing import TestCase, LiveServerTestCase
from pony.orm import db_session, select
from time import sleep
from sources import app, db
import multiprocessing
from database_setup import test_database_create

# multiprocessing.set_start_method("fork")


class AuthTestCase(TestCase):

    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        pass

    def tearDown(self):
        restore_db()

    def register(self, fullname, username, email, password, password2):  # registration function
        return self.client.post(
            '/auth/register',
            data=dict(fullname=fullname, username=username, email=email, password=password, password2=password2,
                      privacy=True, hazard_disclaimer=True),
            follow_redirects=True
        )

    def login(self, username, password):  # login function
        return self.client.post(
            '/auth/login',
            data=dict(username=username, password=password),
            follow_redirects=True
        )

    def password_token_use(self, url): #test the url token
        return self.client.get(
            url,
            follow_redirects=True
        )


    def logout(self):  # logout function
        return self.client.get(
            '/auth/logout',
            follow_redirects=True
        )

    def password_reset_request(self, email):
        return self.client.post(
            '/reset_password_request',
            data=dict(email=email),
            follow_redirects=True
        )

    def email_change_request(self, email, new_email, password):
        return self.client.post(
            '/update_email',
            data=dict(old_password=password, email=email, email2=new_email),
            follow_redirects=True
        )


    def get_emails(self, result_bytes):
        msgs = []  # all the email data are pushed inside an array
        for num in result_bytes[0].split():
            typ, data = con.fetch(num, '(RFC822)')
            msgs.append(data)

        return msgs

    def search(self, key, value, con):
        result, data = con.search(None, key, '"{}"'.format(value))
        return data

    def access_dummy(self, username, password):  # access dummy page function
        """Dummy page is included in the url as I don't know any better
        way of how to make both get and post requests in one method"""

        return self.client.post(
            '/auth/login?next=%2Fhome',
            data=dict(username=username, password=password),
            follow_redirects=True
        )

    def test_successful_user_redirection(self):  # testing successful user redirection to the dummy page
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # registers a new user
        response = self.access_dummy('susan_chem', 'felis catus')  # logs in the registered user
        self.assertIn(b'Welcome to AI4Green, susan_chem!', response.data)  # checks the response message

    def test_successful_user_registration(self):  # testing successful user registration
         response = self.register('tom normal', 'tom_chem', 'tom@mail.com', 'sunny sun', 'sunny sun')  # registers a new user
         self.assertEqual(response.status_code, 200)
         self.assertIn(b'Congratulations, you are now a registered user!', response.data)  # checks the response message

    def test_failed_user_registration_empty_email(self):  # testing failed user registration with empty email field
        response = self.register('tom normal', 'tom_chem', '', 'sunny sun', 'sunny sun')  # registers a new user with the empty email field
        self.assertIn(b'This field is required.', response.data)  # checks the response message

    def test_failed_user_registration_empty_password(self):  # testing failed user registration with empty password
        response = self.register('tom normal', 'tom_chem', 'tom@mail.com', '', 'sunny sun')  # registers a new user with the empty password
        self.assertIn(b'This field is required.', response.data)  # checks the response message

    def test_failed_user_registration_empty_password2(self):  # testing failed user registration with empty password2
        response = self.register('tom normal', 'tom_chem', 'tom@mail.com', 'sunny sun', '')  # registers a new user with the empty password2
        self.assertIn(b'This field is required.', response.data)  # checks the response message

    def test_user_registration_empty_user_type(self):  # testing empty user type as usertype should have a default of 0 and so we should register
        response = self.register('tom normal', 'tom_chem', 'tom@mail.com', 'sunny sun', 'sunny sun')  # registers a new user with the empty password2
        self.assertIn(b'Congratulations, you are now a registered user!', response.data)  # checks the response message

    def test_failed_user_registration_username(self):  # testing failed user registration with a repeating username
        self.register('tom normal', 'tom_chem', 'tom@mail.com', 'sunny sun', 'sunny sun')  # registers a new user
        response = self.register('tom normal', 'tom_chem', 'tom@newmail.com', 'sunny sun', 'sunny sun')  # registers another user with the same name
        self.assertIn(b'Please use a different username.', response.data)  # checks the response message

    def test_failed_user_registration_email(self):  # testing failed user registration with a repeating email
        self.register('tom normal', 'tom_chem', 'tom@mail.com', 'sunny sun', 'sunny sun')  # registers a new user
        response = self.register('tom_normal', 'tim_chem', 'tom@mail.com', 'sunny sun', 'sunny sun')  # registers another user with the same email
        self.assertIn(b'Please use a different email address.', response.data)  # checks the response message

    def test_failed_user_registration(self):  # testing failed user registration with repeating user data
        self.register('tom normal', 'tom_chem', 'tom@mail.com', 'sunny sun', 'sunny sun')  # registers a new user
        response = self.register('tom normal', 'tom_chem', 'tom@mail.com', 'sunny sun', 'sunny sun')  # registers another user with the same data
        self.assertIn(b'Please use a different username.', response.data)  # checks the first response message
        self.assertIn(b'Please use a different email address.', response.data)  # checks the second response message

    def test_successful_user_login(self):  # testing successful user login
        # self.app.get('/auth/register', follow_redirects=True)  # goes to the registration page
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        # self.app.get('/auth/login', follow_redirects=True)  # goes to the login page
        response = self.login('susan_chem', 'felis catus')  # enters the data of the registered user
        #self.assertEqual(response.status_code, 200)  # checks response code
        self.assertIn(b'Welcome to AI4Green, susan_chem!', response.data)  # checks the response message


    def test_failed_user_login_wrong_username(self):  # testing failed user login with wrong username
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        # self.app.get('/auth/login', follow_redirects=True)  # goes to the login page
        response = self.login('sue', 'felis catus')  # enters incorrect username
        # self.assertEqual(response.status_code, 200)  # checks response code
        self.assertIn(b'Invalid username or password', response.data)  # checks the response message

    def test_failed_user_login_wrong_password(self):  # testing failed user login with wrong password
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        response = self.login('susan_chem', 'car')  # enters incorrect password
        # self.assertEqual(response.status_code, 200)  # checks response code
        self.assertIn(b'Invalid username or password', response.data)  # checks the response message

    def test_failed_user_login(self):  # testing failed user login with wrong user data
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        response = self.login('sue', 'car')  # enters incorrect user data
        # self.assertEqual(response.status_code, 200)  # checks response code
        self.assertIn(b'Invalid username or password', response.data)  # checks the response message

    def test_user_logout(self):  # testing user logout
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        self.login('susan_chem', 'felis catus')  # logs in the registered user
        response = self.logout()  # logs out
        self.assertEqual(response.status_code, 200)  # checks the response code
        self.assertIn(b'Please log in', response.data)  # checks the response message

    def test_password_reset_request_redirect(self):  #if you are logged in, reset password redirects to home
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        self.login('susan_chem', 'felis catus')  # logs in the registered user
        response = self.password_reset_request('susan@mail.com')
        self.assertIn(b'Welcome to AI4Green, susan_chem!', response.data)

    def test_failed_registration_short_pw(self):
        """Tests a password must be a certain length"""
        response = self.register('susan smith', 'susan_chem', 'susan@mail.com', 'cat', 'cat')
        self.assertIn(b'Field must be at least 8 characters long.', response.data)

    def test_case_insensitive_username(self):
        """Tests usernames are not case sensitive"""
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        response = self.register('susan smith', 'Susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        self.assertIn(b'Please use a different username.', response.data)

    def test_case_insensitive_email(self):
        """Tests usernames are not case sensitive"""
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        response = self.register('susan smith', 'susan_chem', 'Susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        self.assertIn(b'Please use a different email address.', response.data)

    def test_login_username_case_insensitive(self):
        """Tests usernames are not case sensitive when logging in"""
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        response = self.login('Susan_Chem', 'felis catus')
        self.assertIn(b'Welcome to AI4Green, susan_chem!', response.data)

    def test_username_regex_success(self):
        """Username regex allows only letters, numbers and underscores"""
        response = self.register('susan smith', 'susan_CH3m15ry', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        self.assertIn(b'Congratulations, you are now a registered user!', response.data)

    def test_username_regex_fail(self):
        response = self.register('susan smith', 'sus@NCh£M', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        self.assertIn(b'Invalid input.', response.data)
    #
    # def test_password_reset_request(self):  #if you not logged in, reset password
    #     self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
    #     response_ping = self.password_reset_request('susan@mail.com')
    #
    #     sleep(10)
    #     user = 'AI4Green.alert@gmail.com'
    #     password = 'Ajax2021'
    #     imap_url = 'imap.gmail.com'
    #     con = imaplib.IMAP4_SSL(imap_url)
    #     con.login(user, password)
    #     response, data = con.select('INBOX')
    #     status, data = con.search(None, 'ALL')
    #     msg_list = data[0].split()
    #     msg_list.reverse()
    #     typ, data = con.fetch(msg_list[0], '(RFC822)')
    #     #print('Message %s\n%s\n' % (msg_list[0], data[0][1]))
    #     self.assertIn(b'Failed-Recipients: susan@mail.com', data[0][1])
    #     self.assertIn(b'Check your email for the instructions to reset your password', response_ping.data)
    #     con.store(msg_list[0], '+FLAGS', '\\Deleted')
    #     con.expunge()
    #     con.select('"[Gmail]/Sent Mail"')
    #     status, data = con.search(None, 'ALL')
    #     msg_list = data[0].split()
    #     con.store(msg_list[0], '+FLAGS', '\\Deleted')
    #     con.expunge()
    #     con.close()
    #     con.logout()
    #     sleep(30)
    #
    # def test_password_reset_token(self):  #if you not logged in, reset password
    #     self.register('dan brown', 'dan_chem', 'brown@mail.com', 'felis catus', 'felis catus')  # enters new user data
    #     response_ping = self.password_reset_request('brown@mail.com')
    #
    #     sleep(10)
    #     user = 'AI4Green.alert@gmail.com'
    #     password = 'Ajax2021'
    #     imap_url = 'imap.gmail.com'
    #     con = imaplib.IMAP4_SSL(imap_url, 993)
    #     con.login(user, password)
    #     #print(con.list())
    #
    #     response, data = con.select('"[Gmail]/Sent Mail"')
    #     status, data = con.search(None, 'ALL')
    #     msg_list = data[0].split()
    #
    #     """
    #     The message list is such that the last id is the most recent
    #     """
    #     #msg_list.reverse()
    #     #print(msg_list[0], msg_list[1])
    #     typ, data = con.fetch(msg_list[-1], '(RFC822)')
    #     # print('Message %s\n%s\n' % (msg_list[-1], data[0][1]))
    #     self.assertIn(b'for <brown@mail.com>', data[0][1])
    #     self.assertIn(b'Dear dan_chem', data[0][1])
    #
    #     for response in data:
    #         if isinstance(response, tuple):
    #             msg = email.message_from_bytes(response[1])
    #
    #             subject, encoding = decode_header(msg["Subject"])[0]
    #             if isinstance(subject, bytes):
    #             # if it's a bytes, decode to str
    #                 subject = subject.decode(encoding)
    #             # decode email sender
    #             From, encoding = decode_header(msg.get("From"))[0]
    #             if isinstance(From, bytes):
    #                 From = From.decode(encoding)
    #             # print("Subject:", subject)
    #             # print("From:", From)
    #             # # if the email message is multipart
    #             if msg.is_multipart():
    #             # iterate over email parts
    #                 for part in msg.walk():
    #                     # extract content type of email
    #                     content_type = part.get_content_type()
    #                     content_disposition = str(part.get("Content-Disposition"))
    #                     try:
    #                         # get the email body
    #                         body = part.get_payload(decode=True).decode()
    #                     except:
    #                         pass
    #                     if content_type == "text/plain" and "attachment" not in content_disposition:
    #                         # print text/plain emails and skip attachments
    #                         urls = re.findall('http?://localhost/reset_password/(?:[-\w.]|(?:%[\da-fA-F]{2}))+', body)
    #
    #                     # elif "attachment" in content_disposition:
    #                     #     # download attachment
    #                     #     filename = part.get_filename()
    #                     #     if filename:
    #                     #         folder_name = clean(subject)
    #                     #         if not os.path.isdir(folder_name):
    #                     #         # make a folder for this email (named after the subject)
    #                     #             os.mkdir(folder_name)
    #                     #         filepath = os.path.join(folder_name, filename)
    #                     #         # download attachment and save it
    #                     #         open(filepath, "wb").write(part.get_payload(decode=True))
    #             else:
    #                 # extract content type of email
    #                 content_type = msg.get_content_type()
    #                 # get the email body
    #                 body = msg.get_payload(decode=True).decode()
    #                 if content_type == "text/plain":
    #                     pass
    #                     # print only text email parts
    #                     #print(body)
    #             # if content_type == "text/html":
    #             #     # if it's HTML, create a new HTML file and open it in browser
    #             #     folder_name = clean(subject)
    #             #     if not os.path.isdir(folder_name):
    #             #         # make a folder for this email (named after the subject)
    #             #         os.mkdir(folder_name)
    #             #     filename = "index.html"
    #             #     filepath = os.path.join(folder_name, filename)
    #             #     # write the file
    #             #     open(filepath, "w").write(body)
    #             #     # open in the default browser
    #             #     webbrowser.open(filepath)
    #             # print("=" * 100)
    #
    #     #urls = re.findall('http?://(?:[-\w.]|(?:%[\da-fA-F]{2}))+', msg)
    #
    #     # for line in msg:
    #     #     print(line)
    #     #     urls = re.findall('http?://(?:[-\w.]|(?:%[\da-fA-F]{2}))+', line)
    #     #     print(urls)
    #     # body = msg.get_payload(decode=True).decode
    #     # if content_type == "text/plain":
    #     #     # print only text email parts
    #     #     print(body)
    #     #
    #
    #
    #     #print(urls)
    #     response = self.password_token_use(urls[0])
    #     #print(response.data)
    #     self.assertIn(b'Reset Your Password', response.data)
    #
    #     con.store(msg_list[-1], '+FLAGS', '\\Deleted')
    #     con.expunge()
    #     con.select('INBOX')
    #     status, data = con.search(None, 'ALL')
    #     msg_list = data[0].split()
    #     con.store(msg_list[0], '+FLAGS', '\\Deleted')
    #     con.expunge()
    #     con.close()
    #     con.logout()
    #     sleep(30)


    def test_password_reset_request_fail1(self):  #if you not logged in, reset password with invalid email
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        response = self.password_reset_request('susan@wrongmail')
        self.assertIn(b'Invalid email address', response.data)

    def test_password_reset_request_fail2(self):  #if you not logged in, reset password with invalid email
        self.register('susan smith', 'susan_chem', 'susan@mail.com', 'felis catus', 'felis catus')  # enters new user data
        response = self.password_reset_request('susan')
        self.assertIn(b'Invalid email address', response.data)


# test registration checkboxes
class RegistrationCheckboxesTest(LiveServerTestCase, TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        setup_selenium(self)
        self.driver.find_element_by_id('accept_cookies').click()
        self.driver.find_element_by_id('register').click()

    def tearDown(self):
        self.driver.quit()
        restore_db()

    @db_session
    def test_successful_registration(self):
        self.driver.find_element_by_id('username').send_keys("sunsun")
        self.driver.find_element_by_id('fullname').send_keys("sunsun")
        self.driver.find_element_by_id('email').send_keys("solar@mail.com")
        self.driver.find_element_by_id('password').send_keys("moon moon")
        self.driver.find_element_by_id('password2').send_keys("moon moon")
        self.driver.find_element_by_id('hazard_link').click()
        self.driver.find_element_by_id('privacy_link').click()
        self.driver.find_element_by_id('hazard_disclaimer').click()
        self.driver.find_element_by_id('privacy').click()
        self.driver.find_element_by_id('submit').click()
        sleep(3)
        self.assertIn('Congratulations, you are now a registered user!', self.driver.page_source)
        user = select(x.username for x in db.User if x.email == "solar@mail.com").first()
        self.assertEqual("sunsun", user)

    # def test_hazard_disclaimer_disabled(self):
    #     """Tests the checkbox is initially disabled"""
    #     hazard_disclaimer = self.driver.find_element_by_id('hazard_disclaimer').get_property('disabled')
    #     self.assertEqual(True, hazard_disclaimer)
    #     self.driver.find_element_by_id('hazard_link').click()
    #     self.driver.switch_to.window(self.driver.window_handles[-1])
    #     self.assertIn('<h1>\n        Hazard Disclaimer\n    </h1>', self.driver.page_source)
    #     self.driver.switch_to.window(self.driver.window_handles[0])
    #     hazard_disclaimer = self.driver.find_element_by_id('hazard_disclaimer').get_property('disabled')
    #     self.assertEqual(False, hazard_disclaimer)

    # def test_privacy_disclaimer_disabled(self):
    #     """Tests the checkbox is initially disabled"""
    #     privacy = self.driver.find_element_by_id('privacy').get_property('disabled')
    #     self.assertEqual(True, privacy)
    #     self.driver.find_element_by_id('privacy_link').click()
    #     self.driver.switch_to.window(self.driver.window_handles[-1])
    #     self.assertIn('<h1>\n    AI4Green Privacy Policy\n    </h1>', self.driver.page_source)
    #     self.driver.switch_to.window(self.driver.window_handles[0])
    #     privacy = self.driver.find_element_by_id('privacy').get_property('disabled')
    #     self.assertEqual(False, privacy)

    def test_submit_button(self):
        """Tests submit button is initially disabled"""
        self.driver.find_element_by_id('username').send_keys("sun")
        self.driver.find_element_by_id('fullname').send_keys("sun")
        self.driver.find_element_by_id('email').send_keys("solar@mail.com")
        self.driver.find_element_by_id('password').send_keys("moon")
        self.driver.find_element_by_id('password2').send_keys("moon")
        submit = self.driver.find_element_by_id('submit').get_property('disabled')
        self.assertEqual(True, submit)
        self.driver.find_element_by_id('hazard_link').click()
        self.driver.find_element_by_id('privacy_link').click()
        self.driver.switch_to.window(self.driver.window_handles[0])
        self.driver.find_element_by_id('hazard_disclaimer').click()
        self.driver.find_element_by_id('privacy').click()
        self.driver.find_element_by_id('submit').click()
        submit = self.driver.find_element_by_id('submit').get_property('disabled')
        self.assertEqual(False, submit)


if __name__ == '__main__':
    unittest.main()
