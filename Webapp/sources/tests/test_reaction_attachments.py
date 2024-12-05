from auxiliary_for_tests import *
from sources import app, db
import time
import flask_testing
from database_setup import test_database_create
import unittest
from pathlib import Path
from pony.orm import select, db_session
from datetime import datetime
import pytz


class FileUploadTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
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
        self.current_dir = Path(__file__).parent

    def load_reaction(self):
        select_workgroup(self)
        select_workbook(self, idx=0)
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        time.sleep(8)

    def test_upload_after_delete(self):
        with db_session:
            reaction = select(x for x in db.Reaction if x.name == 'nitro reduction2').first()
            current_time = datetime.now(pytz.timezone('Europe/London')).replace(tzinfo=None)
            file_details = {'dummy': 'dummy'}
            with db_session:
                for x in range(5):
                    db.ReactionDataFile(reaction=reaction, storage_name='dummy', container_name='dummy',
                                        uuid=f'dummy{x}', display_name='dummy', time_of_upload=current_time,
                                        file_details=file_details)
            time.sleep(1)
        self.load_reaction()
        for i in range(5):
            try:
                self.driver.find_element_by_id('delete-file-button1').click()
                time.sleep(2)
                self.driver.switch_to.alert.accept()
                time.sleep(1)
                break
            except:
                pass
        test_image = str(self.current_dir / 'test_data' / 'nmr.png')
        self.upload_file(test_image)
        attachment_name = self.driver.find_element_by_id('view-file-button5').text
        self.assertEqual('nmr.png', attachment_name)


    def test_attachment_appears(self):
        self.load_reaction()
        test_image = str(self.current_dir / 'test_data' / 'nmr.png')
        self.upload_file(test_image)
        attachment = self.driver.find_element_by_id('view-file-button1')
        self.assertIsNotNone(attachment)

    def upload_file(self, file):
        web_upload = self.driver.find_element_by_id('upload-files')
        web_upload.send_keys(file)
        time.sleep(1)
        for i in range(5):
            try:
                self.driver.find_element_by_id('file-upload').click()
                time.sleep(3)
                self.driver.switch_to.alert.dismiss()
                time.sleep(2)
                break
            except:
                pass


    def test_mimetype(self):
        """Text.png is actually a .txt file - should be rejected by backend"""
        self.load_reaction()
        test_image = str(self.current_dir / 'test_data' / 'text.png')
        web_upload = self.driver.find_element_by_id('upload-files')
        web_upload.send_keys(test_image)
        time.sleep(2)
        self.driver.find_element_by_id('file-upload').click()
        time.sleep(2)
        try:
            self.driver.find_element_by_id('view-file-button1')
            not_found = False
        except:
            not_found = True
        self.assertTrue(not_found)

    def test_max_files(self):
        with db_session:
            reaction = select(x for x in db.Reaction if x.name == 'nitro reduction2').first()
            current_time = datetime.now(pytz.timezone('Europe/London')).replace(tzinfo=None)
            file_details = {'dummy': 'dummy'}
            with db_session:
                for x in range(10):
                    db.ReactionDataFile(reaction=reaction, storage_name='dummy', container_name='dummy',
                                        uuid=f'dummy{x}', display_name='dummy', time_of_upload=current_time,
                                        file_details=file_details)
        time.sleep(1)
        self.load_reaction()
        # reaction has 10, no more should be uploaded
        test_image = str(self.current_dir / 'test_data' / 'nmr.png')
        web_upload = self.driver.find_element_by_id('upload-files')
        web_upload.send_keys(test_image)
        for i in range(5):
            try:
                self.driver.find_element_by_id('file-upload').click()
                time.sleep(2)
                alert_text = self.driver.switch_to.alert.text
                break
            except:
                pass
        self.assertEqual("The maximum number of attachments for a reaction is 10.", alert_text)

    def test_size_limit(self):
        self.load_reaction()
        test_large_file = str(self.current_dir / 'test_data' / 'selenium_tutorial.pdf')
        self.upload_file(test_large_file)
        try:
            self.driver.find_element_by_id('view-file-button1')
            not_found = False
        except:
            not_found = True
        self.assertTrue(not_found)

    def test_pdf(self):
        self.load_reaction()
        test_pdf = str(self.current_dir / 'test_data' / 'AI4Green_quick_guide.pdf')
        self.upload_file(test_pdf)
        attachment = self.driver.find_element_by_id('view-file-button1')
        self.assertIsNotNone(attachment)

    def test_file_delete(self):
        self.load_reaction()
        test_image = str(self.current_dir / 'test_data' / 'nmr.png')
        self.upload_file(test_image)
        attachment = self.driver.find_element_by_id('view-file-button1')
        self.assertIsNotNone(attachment)
        self.driver.find_element_by_id('delete-file-button1').click()
        time.sleep(2)
        self.driver.switch_to.alert.accept()
        time.sleep(1)
        try:
            self.driver.find_element_by_id('view-file-button1')
            not_found = False
        except:
            not_found = True
        self.assertTrue(not_found)


    def test_file_view(self):
        self.load_reaction()
        test_image = str(self.current_dir / 'test_data' / 'nmr.png')
        self.upload_file(test_image)
        self.driver.find_element_by_id('view-file-button1').click()
        time.sleep(2)
        self.driver.switch_to.window(self.driver.window_handles[1])
        self.assertIn('data:image/png;base64,iVBORw0KGgoAAAANS' ,self.driver.page_source)


class FileDownloadTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def setUp(self):
        """We load a test instance of the app, clear the db, register a new user then login"""
        self.current_dir = Path(__file__).parent
        options = webdriver.ChromeOptions()
        options.add_experimental_option("prefs", {
            "download.default_directory": str(self.current_dir / 'test_data' / 'test_downloads'),
            "download.prompt_for_download": False,
            "download.directory_upgrade": True,
            "safebrowsing.enabled": True
        })
        chrome_dir = os.path.join(basedir, "sources/tests/test_data/chromedriver")
        options.add_argument("--start-maximized")
        self.driver = webdriver.Chrome(chrome_dir, options=options)
        self.driver.get("http://localhost:8943")
        login(self, "BB_Test", "BB_login")

    def load_reaction(self):
        select_workgroup(self)
        select_workbook(self, idx=0)
        self.driver.find_element_by_id('nitro reduction2-reload').click()
        time.sleep(8)

    def test_download_file(self):
        self.load_reaction()
        test_image = str(self.current_dir / 'test_data' / 'nmr.png')
        web_upload = self.driver.find_element_by_id('upload-files')
        web_upload.send_keys(test_image)
        time.sleep(2)
        for i in range(5):
            try:
                self.driver.find_element_by_id('file-upload').click()
                time.sleep(3)
                self.driver.switch_to.alert.dismiss()
                time.sleep(2)
                break
            except:
                pass
        for i in range(5):
            try:
                self.driver.find_element_by_id('download-file-button1').click()
                time.sleep(2)
                file = self.current_dir / 'test_data' / 'test_downloads' / 'nmr.png'
                x = Path.exists(file)
                self.assertTrue(x)
                break
            except:
                pass
        file.unlink()
