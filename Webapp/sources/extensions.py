"""
Extensions module.

"""

from flask_login import LoginManager
from flask_marshmallow import Marshmallow
from flask_migrate import Migrate
from flask_sqlalchemy import SQLAlchemy
from sources.email_sender import EmailSender

login = LoginManager()
mail = EmailSender()
ma = Marshmallow()
db = SQLAlchemy()
migrate = Migrate()
