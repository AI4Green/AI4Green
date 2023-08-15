from auxiliary_for_tests import *
import unittest
from selenium import webdriver
import time
import flask_testing
from selenium.webdriver.common.keys import Keys

from utilities import basedir
from sources import app, db
import os
from database_setup import test_database_create


# test news feed
class NewsFeedTest(flask_testing.LiveServerTestCase, flask_testing.TestCase):
    def create_app(self):
        app.config.from_object('config.TestConfig')
        test_database_create(db)
        return app

    def setUp(self):
        self.login()

    def tearDown(self):
        self.driver.quit()
        restore_db()

    def login(self):
        """This function sets the headless browser up, logs the user in and is called from setUp"""
        setup_selenium(self)
        login(self, "admin", "admin_login")

    def post_message(self, title, message):
        self.driver.find_element_by_id('news-feed-tab').click()
        time.sleep(1)
        news_post_title = self.driver.find_element_by_id('newsPostTitle')
        news_post_title.clear()
        news_post_title.send_keys(title)
        news_post_message = self.driver.find_element_by_id('newsPostMessage')
        news_post_message.clear()
        news_post_message.send_keys(message)
        self.driver.find_element_by_id('submit-post').click()
        # test success message
        self.assertIn("Your message has been posted!", self.driver.page_source)

    # test news feed
    def test_news_feed(self):
        # test initial message there are no news items
        self.assertIn("No news items to show!", self.driver.page_source)
        # add two news item
        self.driver.find_element_by_id('admin_dashboard').click()
        self.post_message("First news item!", "Some content.")
        self.post_message("Second news item!", "Some different content.")
        self.driver.find_element_by_id("TopNavHomeButton").click()
        # test initial message is not there
        self.assertNotIn("No news items to show!", self.driver.page_source)
        # test new items appear with most recent first
        self.assertIn("First news item", self.driver.page_source)
        self.assertTrue(self.driver.page_source.find('Second news item') <
                        self.driver.page_source.find('First news item'))
        # test can delete first
        self.driver.find_element_by_id("delete-1").click()
        self.assertIn("This message has been deleted!", self.driver.page_source)
        self.assertNotIn("First news item", self.driver.page_source)
        self.assertIn("delete-2", self.driver.page_source)
        # test non-admin cannot delete
        logout(self)
        login(self, "SR_Test", "SR_login")
        self.assertNotIn("delete-2", self.driver.page_source)

