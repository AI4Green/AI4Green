from flask import Blueprint

news_feed_bp = Blueprint("news_feed", __name__)

from . import routes
