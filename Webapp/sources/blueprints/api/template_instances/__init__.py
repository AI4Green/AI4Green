from apiflask import APIBlueprint

template_instances_api_bp = APIBlueprint(
    "template_instances", __name__, url_prefix="/template-instances"
)

from . import routes
