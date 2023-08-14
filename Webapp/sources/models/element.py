from sources.extensions import db
from .base import Model


class Element(Model):
    __tablename__ = "Element"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    name = db.Column(db.Text, nullable=False)
    symbol = db.Column(db.Text, nullable=False)
    remaining_supply = db.Column(db.Text, nullable=False)
    colour = db.Column(db.Text, nullable=False)
