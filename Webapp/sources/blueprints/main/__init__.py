# from flask import Blueprint
from apiflask import APIBlueprint

main_bp = APIBlueprint("main", __name__)

from . import routes
