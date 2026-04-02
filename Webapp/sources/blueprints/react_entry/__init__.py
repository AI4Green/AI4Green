from apiflask import APIBlueprint

# spa = single page application
# this route ports the front end to react instead of flask/jinja
spa_bp = APIBlueprint("spa", __name__, url_prefix="/spa")

from . import routes
