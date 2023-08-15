from sources.extensions import db
from .base import Model


class NewsItem(Model):
    __tablename__ = "NewsItem"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    title = db.Column(db.Text, nullable=False)
    message = db.Column(db.Text, nullable=False)
    time = db.Column(db.DateTime, nullable=False)
