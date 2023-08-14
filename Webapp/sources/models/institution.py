from sources.extensions import db
from .base import Model


class Institution(Model):
    __tablename__ = "Institution"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    name = db.Column(db.Text, nullable=False, unique=True)
    time_of_creation = db.Column(db.DateTime)
