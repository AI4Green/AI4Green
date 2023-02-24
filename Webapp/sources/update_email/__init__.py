from flask import Blueprint

update_email_bp = Blueprint('update_email', __name__)

from sources.update_email import routes
