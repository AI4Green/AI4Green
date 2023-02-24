#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The main file of the app
"""

from flask import Flask  # imports an instance of Flask
from pony.orm import Database  # Imports the class Database
from pony.flask import Pony  # Imports Pony which manages the db_sessions for view functions
from flask_login import LoginManager  # imports the login manager from Flask-Login
from flask_mail import Mail  # import mail service from flask
# The extension Flask-Login manages the user logged-in state
#import config  # imports the Config class
import os


def create_application(is_test):
    """The function takes the variable from the environmental variable. The value 1 means the app will start in testing/development
    the value 0 or no value means the app will start in deployment mode."""
    app = Flask(__name__)
    if is_test == '0':
        db_name = 'PONY_DATABASE'
        print("Starting app with postgres database")
        app.config.from_object('config.DeploymentConfig')  # tells Flask to read config.py and apply it
    elif is_test == '1':
        db_name = 'PONY_TEST'
        print("Starting app in test/development mode with sqlite database")
        app.config.from_object('config.TestConfig')  # tells Flask to read config.py and apply it
    else:
        raise Exception("Invalid value assigned to the 'TESTING' environmental variable. Assigned value must be '1' or '0'")
    return app, db_name


# get the value for the environmental variable
is_test_env = os.environ.get("TESTING")
# if the environmental variable hasn't been assigned, it will default to 1 for testing/development mode
is_test = '1' if is_test_env is None else is_test_env
# an instance of the app is made and the app config including is assigned values. db_name specifies the database to bind.
app, db_name = create_application(is_test)
print(db_name)
db = Database()  # Creates an instance of a Pony database
Pony(app)  # links the app to the user database
login = LoginManager(app)  # initializes the login manager
login.login_view = 'auth.login'  # the view function that handles logins
# This function forces users to log in before they can view certain pages of the app

# set up mail service
mail = Mail(app)

from sources import models  # Importing the models module - contains all database tables
db.bind(**app.config[db_name])  # Binds the database using the config app
models.define_database(db)
db.generate_mapping(create_tables=True)  # Makes tables based on the structure imported from models

"""Now we need to import the app modules as blueprints"""
from sources.workgroup import workgroup_bp
app.register_blueprint(workgroup_bp)

from sources.main import main_bp
app.register_blueprint(main_bp)

from sources.auth import auth_bp
app.register_blueprint(auth_bp, url_prefix='/auth')

from sources.reaction_table import reaction_table_bp
app.register_blueprint(reaction_table_bp)

from sources.reagents import reagents_bp
app.register_blueprint(reagents_bp)

from sources.notifications import notifications_bp
app.register_blueprint(notifications_bp)

from sources.solvents import solvents_bp
app.register_blueprint(solvents_bp)

from sources.summary import summary_bp
app.register_blueprint(summary_bp)

from sources.standard_landing_page import standard_landing_page_bp
app.register_blueprint(standard_landing_page_bp)

from sources.delete_profile import delete_profile_bp
app.register_blueprint(delete_profile_bp)

from sources.join_workgroup import join_workgroup_bp
app.register_blueprint(join_workgroup_bp)

from sources.create_workgroup import create_workgroup_bp
app.register_blueprint(create_workgroup_bp)

from sources.create_workbook import create_workbook_bp
app.register_blueprint(create_workbook_bp)

from sources.manage_workgroup import manage_workgroup_bp
app.register_blueprint(manage_workgroup_bp)

from sources.manage_workbook import manage_workbook_bp
app.register_blueprint(manage_workbook_bp)

from sources.update_email import update_email_bp
app.register_blueprint(update_email_bp)

from sources.workgroup_membership_summary import workgroup_membership_summary_bp
app.register_blueprint(workgroup_membership_summary_bp)

from sources.save_reaction import save_reaction_bp
app.register_blueprint(save_reaction_bp)

from sources.reset_password import reset_password_bp
app.register_blueprint(reset_password_bp)

from sources.novel_compound import novel_compound_bp
app.register_blueprint(novel_compound_bp)

from sources.compound_data_error_report import compound_data_error_report_bp
app.register_blueprint(compound_data_error_report_bp)

from sources.admin_dashboard import admin_dashboard_bp
app.register_blueprint(admin_dashboard_bp)

from sources.solvent_guide import solvent_guide_bp
app.register_blueprint(solvent_guide_bp)

from sources.news_feed import news_feed_bp
app.register_blueprint(news_feed_bp)

from sources.export_data import export_data_bp
app.register_blueprint(export_data_bp)
""""""

app.config["PROPAGATE_EXCEPTIONS"] = False


@app.context_processor
def inject_session_context():
    """This hides the footer in testing mode through a jinja variable sent to base html for every rendered template"""
    if 'UNIT_TEST' in os.environ:
        return dict(session_type="UNIT_TEST", marvin_js_key=app.config["MARVIN_JS_API_KEY"])
    else:
        return dict(session_type="", marvin_js_key=app.config["MARVIN_JS_API_KEY"])


if __name__ == "__main__":
    app.run()  # to run the app in PyCharm
