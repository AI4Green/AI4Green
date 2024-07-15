from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model

metadata = db.Model.metadata


class Retrosynthesis(Model):
    def __init__(self, **kwargs) -> None:
        self.time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        super().__init__(**kwargs)

    __tablename__ = "Retrosynthesis"
    __table_args__ = (db.UniqueConstraint("workbook", "name"),)

    id = db.Column(db.Integer, primary_key=True)

    name = db.Column(db.Text, nullable=False)
    target_smiles = db.Column(db.Text, nullable=False)
    creator = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    creator_person = db.relationship("Person", foreign_keys=[creator])

    time_of_creation = db.Column(db.DateTime, nullable=False)

    routes = db.Column(db.JSON, nullable=False)
    conditions = db.Column(db.JSON, nullable=True)
    sustainability = db.Column(db.JSON, nullable=True)
    uuid = db.Column(db.Text, nullable=False)

    workbook = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False
    )
    Workbook = db.relationship("WorkBook")
