from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model


class PolymerRepeatUnit(Model):
    __tablename__ = "PolymerRepeatUnit"

    def __init__(self, **kwargs) -> None:
        self.time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    polymer_id = db.Column(
        db.ForeignKey("PolymerNovelCompound.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    polymer = db.relationship("PolymerNovelCompound", backref="PolymerRepeatUnit")
    smiles = db.Column(db.Text, nullable=False)
    molec_weight = db.Column(db.Float(53))
    molec_formula = db.Column(db.Text, nullable=False)
    workbook = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False, index=True
    )
    time_of_creation = db.Column(db.DateTime)
