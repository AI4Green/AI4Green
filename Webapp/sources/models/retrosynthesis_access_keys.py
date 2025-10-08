from uuid import uuid4

from sources.extensions import db
from sqlalchemy.types import Uuid

from .base import Model


class RetrosynthesisKey(Model):
    __tablename__ = "RetrosynthesisKey"

    id = db.Column(db.Integer, primary_key=True)
    key = db.Column(Uuid, default=uuid4, nullable=False, unique=True)

    # one-to-one relationship to user
    user_id = db.Column(db.ForeignKey("User.id"), nullable=True, unique=True)
    user = db.relationship("User", back_populates="RetrosynthesisKey")
