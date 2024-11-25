#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The main file of the app
"""
import os
from typing import Dict

from flask import Flask
from sources import config, models
from sources.auxiliary import get_notification_number, get_workgroups
from sources.extensions import db, login, ma, mail, migrate


def create_app(c: str = "dev") -> Flask:
    """
    Create the Flask app.

    Args:
        config: prod | dev Configuration to use.
        prod means the app will start in deployment mode
        dev or none means the app will start in testing/development

    Returns:
        Flask app.

    """
    app = Flask(__name__)

    if c == "prod":
        print("Starting AI4Green in deployment configuration")
        app.config.from_object(config.DeploymentConfig)
    elif c == "dev":
        print("Starting AI4Green in development configuration")
        app.config.from_object(config.DevConfig)
    else:
        print("Starting AI4Green in testing configuration")
        app.config.from_object(config.TestConfig)

    register_extensions(app)

    app.context_processor(lambda: inject_session_context(app))

    with app.app_context():
        from sources.services.retrosynthesis.dashboard import init_dashboard

        app = init_dashboard(app)
        return run_app(app)


def run_app(app: Flask) -> Flask:
    """
    Create the Flask app with app context.

    Args:
        app: Flask app.

    """
    register_cli(app)
    register_blueprints(app)

    return app


def register_extensions(app: Flask) -> None:
    """
    Registers extensions for the app.

    Args:
        app: Flask app.

    """
    db.init_app(app)
    ma.init_app(app)
    migrate.init_app(app, db)

    login.init_app(app)
    login.login_view = "auth.login"

    @login.user_loader
    def load_user(user_id: int):
        return models.User.query.get(user_id)

    mail.init_app(app)

    return None


def register_cli(app: Flask) -> None:
    """
    Register Flask CLI commands.

    Args:
        app: Flask app.

    """
    from sources.commands import commands

    app.cli.add_command(commands.drop_db)
    app.cli.add_command(commands.create_db)
    app.cli.add_command(commands.seed_db)
    app.cli.add_command(commands.seed_users)
    app.cli.add_command(commands.download_pubchem)
    app.cli.add_command(commands.update_db)
    app.cli.add_command(commands.add_user)


def register_blueprints(app: Flask) -> None:
    """
    Registers blueprints for the app.
    """

    from sources.blueprints.workgroup import workgroup_bp

    app.register_blueprint(workgroup_bp)

    from sources.blueprints.main import main_bp

    app.register_blueprint(main_bp)

    from sources.blueprints.auth import auth_bp

    app.register_blueprint(auth_bp, url_prefix="/auth")

    from sources.blueprints.reaction_table import reaction_table_bp

    app.register_blueprint(reaction_table_bp)

    from sources.blueprints.reagents import reagents_bp

    app.register_blueprint(reagents_bp)

    from sources.blueprints.notifications import notifications_bp

    app.register_blueprint(notifications_bp)

    from sources.blueprints.solvents import solvents_bp

    app.register_blueprint(solvents_bp)

    from sources.blueprints.solvent_PCA import solvent_PCA_bp

    app.register_blueprint(solvent_PCA_bp)

    from sources.blueprints.summary import summary_bp

    app.register_blueprint(summary_bp)

    from sources.blueprints.reaction_list import reaction_list_bp

    app.register_blueprint(reaction_list_bp)

    from sources.blueprints.delete_profile import delete_profile_bp

    app.register_blueprint(delete_profile_bp)

    from sources.blueprints.join_workgroup import join_workgroup_bp

    app.register_blueprint(join_workgroup_bp)

    from sources.blueprints.create_workgroup import create_workgroup_bp

    app.register_blueprint(create_workgroup_bp)

    from sources.blueprints.create_workbook import create_workbook_bp

    app.register_blueprint(create_workbook_bp)

    from sources.blueprints.manage_workgroup import manage_workgroup_bp

    app.register_blueprint(manage_workgroup_bp)

    from sources.blueprints.manage_workbook import manage_workbook_bp

    app.register_blueprint(manage_workbook_bp)

    from sources.blueprints.update_email import update_email_bp

    app.register_blueprint(update_email_bp)

    from sources.blueprints.workgroup_membership_summary import (
        workgroup_membership_summary_bp,
    )

    app.register_blueprint(workgroup_membership_summary_bp)

    from sources.blueprints.save_reaction import save_reaction_bp

    app.register_blueprint(save_reaction_bp)

    from sources.blueprints.reset_password import reset_password_bp

    app.register_blueprint(reset_password_bp)

    from sources.blueprints.email_verification import email_verification_bp

    app.register_blueprint(email_verification_bp)

    from sources.blueprints.novel_compound import novel_compound_bp

    app.register_blueprint(novel_compound_bp)

    from sources.blueprints.compound_data_error_report import (
        compound_data_error_report_bp,
    )

    app.register_blueprint(compound_data_error_report_bp)

    from sources.blueprints.admin_dashboard import admin_dashboard_bp

    app.register_blueprint(admin_dashboard_bp)

    from sources.blueprints.solvent_guide import solvent_guide_bp

    app.register_blueprint(solvent_guide_bp)

    from sources.blueprints.news_feed import news_feed_bp

    app.register_blueprint(news_feed_bp)

    from sources.blueprints.export_data import export_data_bp

    app.register_blueprint(export_data_bp)

    from sources.blueprints.search import search_bp

    app.register_blueprint(search_bp)

    from sources.blueprints.retrosynthesis import retrosynthesis_bp

    app.register_blueprint(retrosynthesis_bp)

    from sources.blueprints.version import version_bp

    app.register_blueprint(version_bp)

    from sources.blueprints.utils import utils_bp

    app.register_blueprint(utils_bp)


def inject_session_context(app: Flask) -> Dict[str, str]:
    """
    This hides the footer in testing mode through a jinja variable sent to base html for every rendered template

    """
    workgroups = get_workgroups()
    notification_number = get_notification_number()
    if "UNIT_TEST" in os.environ:
        return dict(
            session_type="UNIT_TEST",
            marvin_js_key=app.config["MARVIN_JS_API_KEY"],
            workgroups=workgroups,
            notification_number=notification_number,
        )
    else:
        return dict(
            session_type="",
            marvin_js_key=app.config["MARVIN_JS_API_KEY"],
            workgroups=workgroups,
            notification_number=notification_number
        )
