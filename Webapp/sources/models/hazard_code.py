from sources.extensions import db
from .base import Model


class HazardCode(Model):
    __tablename__ = "HazardCode"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    code = db.Column(db.Text, nullable=False)
    phrase = db.Column(db.Text, nullable=False)
    category = db.Column(db.Text)
