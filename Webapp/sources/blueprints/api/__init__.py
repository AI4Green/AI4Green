from apiflask import APIBlueprint

from .fields import fields_api_bp
from .input_types import input_types_api_bp
from .sections import sections_api_bp
from .template_instances import template_instances_api_bp
from .templates import templates_api_bp
from .users import users_api_bp

# Define the Master API Blueprint
api = APIBlueprint("api_v1", __name__, url_prefix="/api")

api.register_blueprint(fields_api_bp)
api.register_blueprint(input_types_api_bp)
api.register_blueprint(templates_api_bp)
api.register_blueprint(template_instances_api_bp)
api.register_blueprint(sections_api_bp)
api.register_blueprint(users_api_bp)
