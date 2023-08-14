from sources.extensions import db
from .base import Model


class ReactionNote(Model):
    __tablename__ = "ReactionNote"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    time_of_creation = db.Column(db.DateTime, nullable=False)
    text = db.Column(db.Text, nullable=False)
    author = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    reaction = db.Column(
        db.ForeignKey("Reaction.id", ondelete="CASCADE"), nullable=False, index=True
    )
    Author = db.relationship("Person", backref="ReactionNote")
