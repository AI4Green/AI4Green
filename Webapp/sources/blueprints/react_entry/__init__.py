from apiflask import APIBlueprint

# spa = single page application
# this route ports the front end to react instead of flask/jinja
react_entry_bp = APIBlueprint("react_entry", __name__, url_prefix="/spa")

from . import routes
