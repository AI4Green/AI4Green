#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a configuration module for the app
"""
import os
from utilities import read_yaml, read_yaml_db_config


class BaseConfig(object):  # class to store configuration variables
    """!!!Change the secret key when deploying the app on a server!!!"""
    """Flask and some of its extensions use the value of the secret
    key as a cryptographic key to generate signatures or tokens.
    The Flask-WTF extension uses it to protect web forms against
    Cross-Site Request Forgery."""
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'
    """The Pony from pony.flask extension takes the location of the application's 
    database from the PONY configuration variable.
    The database URL is taken from the DATABASE_URL environment variable, 
    and if that isn't defined, it is connected to a database named users.db 
    located in the main directory of the application, which is stored in the basedir variable."""

    """"""
    DEBUG = False
    WTF_CSRF_ENABLED = False
    LIVESERVER_PORT = 8943
    LIVESERVER_TIMEOUT = 10
    """
    Settings for the mail app
    """
    MAIL_SERVER = read_yaml(["mail", "mail_server"])
    MAIL_PORT = read_yaml(["mail", "mail_port"])
    MAIL_USERNAME = read_yaml(["mail", "mail_username"])
    MAIL_PASSWORD = read_yaml(["mail", "mail_password"])
    MAIL_USE_TLS = read_yaml(["mail", "mail_use_tls"])
    MAIL_USE_SSL = read_yaml(["mail", "mail_use_ssl"])
    MARVIN_JS_API_KEY = read_yaml(['marvin_js_api_key'])


class DeploymentConfig(BaseConfig):
    TESTING = False
    PONY_DATABASE = read_yaml_db_config("db_postgres")


class TestConfig(BaseConfig):
    TESTING = True  #####change to false to get emails to send
    LOGIN_DISABLED = False
    PONY_TEST = read_yaml_db_config("db_sqlite")

