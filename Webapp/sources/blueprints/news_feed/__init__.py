from apiflask import APIBlueprint

news_feed_bp = APIBlueprint("news_feed", __name__)

from . import routes
