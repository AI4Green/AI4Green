from flask import Blueprint

news_feed_bp = Blueprint('news_feed', __name__)

from sources.news_feed import routes
