#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a configuration module for the app
"""
import os

from dotenv import load_dotenv

load_dotenv()


class BaseConfig(object):  # class to store configuration variables
    """!!!Change the secret key when deploying the app on a server!!!"""

    """Flask and some of its extensions use the value of the secret
    key as a cryptographic key to generate signatures or tokens.
    The Flask-WTF extension uses it to protect web forms against
    Cross-Site Request Forgery."""
    SERVER_NAME = os.getenv("SERVER_NAME", "127.0.0.1:80")
    SECRET_KEY = os.getenv("SECRET_KEY", "change-me")
    WTF_CSRF_ENABLED = False
    LIVESERVER_TIMEOUT = 10
    MAX_CONTENT_LENGTH = 1024 * 1024
    UPLOAD_EXTENSIONS = [
        ".arw", "cif", ".cdf", ".cdx", ".csv", ".D", ".dat", ".DTA", ".dx", ".fid", ".gz", ".jpg", ".jdx", ".jcamp",
        ".mnova", ".pdf",  ".pkl", ".png", ".pdb", ".pssession", ".qgd", ".raw", ".res", ".topspin", ".txt", ".vrml",
        ".xlsx", "xls",  ".xyz", ".zip"
    ]
    UPLOAD_MIME_TYPES = [
                        "image/x-sony-arw",
                        "chemical/x-cif",
                        "application/x-netcdf",
                        "chemical/x-cdx",
                        "text/csv",
                        "text/plain",
                        "application/octet-stream",
                        "chemical/x-jcamp-dx",
                        "application/gzip",
                        "image/jpeg",
                        "application/pdf",
                        "image/png",
                        "chemical/x-pdb",
                        "model/vrml",
                        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        "application/vnd.ms-excel",
                        "chemical/x-xyz",
                        "application/zip"
                        ]


    # the number of molecules to take from PubChem
    COMPOUND_LIMIT = os.getenv("COMPOUND_LIMIT", 10000)

    # Mail app - use example values in comments to help setup your own mail server if desired
    MAIL_SERVER = os.getenv("MAIL_SERVER", "None")  # 'smtp.gmail.com'
    MAIL_PORT = os.getenv("MAIL_PORT", "None")  # 465
    MAIL_USERNAME = os.getenv("MAIL_USERNAME", "None")  # mail@.com
    MAIL_PASSWORD = os.getenv("MAIL_PASSWORD", "None")
    MAIL_USE_TLS = os.getenv("MAIL_USE_TLS", False)
    MAIL_USE_SSL = os.getenv("MAIL_USE_SSL", True)
    MAIL_USE_LOCAL = os.getenv("MAIL_USE_LOCAL", "local")
    MAIL_ADMIN_SENDER = os.getenv("MAIL_ADMIN_SENDER", "admin@ai4green.app")
    MAIL_SAVE_DIR = os.getenv("MAIL_SAVE_DIR", "temp")

    # Marvin JS - Change to your own
    MARVIN_JS_API_KEY = os.getenv("MARVIN_JS_API_KEY", "")

    # reCaptcha Keys - Change to your own
    RECAPTCHA_PUBLIC_KEY = os.getenv("RECAPTCHA_PUBLIC_KEY", "")
    RECAPTCHA_PRIVATE_KEY = os.getenv("RECAPTCHA_PRIVATE_KEY", "")

    # Azurite file storage
    AZURE_STORAGE_CONNECTION_STRING = os.getenv(
        "AZURE_STORAGE_CONNECTION_STRING",
        "DefaultEndpointsProtocol=http;AccountName=devstoreaccount1;AccountKey=Eby8vdM02xNOcqFlqUwJPLlmEtlCDXJ1OUzFT50uSRZ6IFsuFq2UVErCz4I6tq/K1SZFPTOtr/KBHBeksoGMGw==;BlobEndpoint=http://127.0.0.1:10000/devstoreaccount1;",
    )
    STORAGE_CONTAINER = os.getenv("STORAGE_CONTAINER", "experiment-data")

    SQLALCHEMY_DATABASE_URI = os.getenv(
        "DATABASE_URL", "postgresql://postgres:postgres@localhost:5432/ai4green"
    )
    SQLALCHEMY_BINDS = {
        "db": SQLALCHEMY_DATABASE_URI,
        "update": "sqlite:///temp_update.sqlite",
    }

    APP_DIRECTORY = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", ".."))

    # Traverse up directories to reach the project root
    PROJECT_ROOT = os.path.abspath(os.path.join(APP_DIRECTORY, ".."))

    # Construct paths relative to the project root
    PYPROJECT_PATH = os.path.join(PROJECT_ROOT, "pyproject.toml")
    HASH_FILE_PATH = os.path.join(PROJECT_ROOT, "hash.txt")

    RETROSYNTHESIS_API_URL = os.getenv(
        "RETROSYNTHESIS_API_URL", "http://127.0.0.1:8000/"
    )
    RETROSYNTHESIS_API_KEY = os.getenv("RETROSYNTHESIS_API_KEY", "retro_key")

    CONDITIONS_API_URL = os.getenv("CONDITIONS_API_URL", "http://127.0.0.1:9901")


class TestConfig(BaseConfig):
    TESTING = True
    DEBUG = True
    LIVESERVER_PORT = 8943
    LOGIN_DISABLED = False
    SQLALCHEMY_DATABASE_URI = os.getenv(
        "TEST_DATABASE_URL",
        "postgresql://postgres:postgres@localhost:5433/ai4greentest",
    )
    SQLALCHEMY_BINDS = {
        "db": SQLALCHEMY_DATABASE_URI,
        "update": "sqlite:///temp_update.sqlite",
    }


class DevConfig(BaseConfig):
    TESTING = False
    DEBUG = True


class DeploymentConfig(BaseConfig):
    TESTING = False
    DEBUG = False
    AZURE_BLOB_SERVICE_URL = os.getenv(
        "AZURE_BLOB_SERVICE_URL", "ai4greenreservestorage"
    )
