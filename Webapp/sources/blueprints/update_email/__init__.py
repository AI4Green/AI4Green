from apiflask import APIBlueprint

update_email_bp = APIBlueprint("update_email", __name__)

from . import routes
