from apiflask import APIBlueprint

email_verification_bp = APIBlueprint("email_verification", __name__)

from . import routes
