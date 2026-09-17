from apiflask import APIBlueprint

from .coshh import coshh_api_bp
from .users import users_api_bp

# Define the Master API Blueprint
api = APIBlueprint("api_v1", __name__, url_prefix="/api")

api.register_blueprint(coshh_api_bp)
api.register_blueprint(users_api_bp)
