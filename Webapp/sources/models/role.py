from sources.extensions import db
from .base import Model


class Role(Model):
    __tablename__ = "Role"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    name = db.Column(db.Text, nullable=False)
    role_description = db.Column(db.Text, nullable=False)
