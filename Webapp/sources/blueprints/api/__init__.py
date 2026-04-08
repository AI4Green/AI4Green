from apiflask import APIBlueprint

from .coshh import coshh_api_bp

# Define the Master API Blueprint
api_v1 = APIBlueprint("api_v1", __name__, url_prefix="/api")

api_v1.register_blueprint(coshh_api_bp)
