from apiflask import APIBlueprint

from .sections import sections_api_bp
from .templates import templates_api_bp
from .users import users_api_bp

# Define the Master API Blueprint
api = APIBlueprint("api_v1", __name__, url_prefix="/api")

api.register_blueprint(templates_api_bp)
api.register_blueprint(sections_api_bp)
api.register_blueprint(users_api_bp)
