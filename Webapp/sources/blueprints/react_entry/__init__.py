from apiflask import APIBlueprint

# spa = single page application
react_entry_bp = APIBlueprint("react_entry", __name__, url_prefix="/spa")

from . import routes
