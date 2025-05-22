from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model


class ReactionEditingHistory(Model):
    __tablename__ = "reaction_editing_history"

    def __init__(self, **kwargs) -> None:
        self.time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)
    workbook_id = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=True, index=True
    )
    reaction_id = db.Column(db.Integer, nullable=False, index=True)
    person_id = db.Column(db.Integer, db.ForeignKey("Person.id"), nullable=False)
    time = db.Column(db.DateTime, nullable=False)
    field_name = db.Column(db.Text, nullable=False)
    change_details = db.Column(db.JSON)

    person = db.relationship("Person", backref="ReactionEditingHistory")
