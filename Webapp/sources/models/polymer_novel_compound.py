from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model


class PolymerNovelCompound(Model):
    __tablename__ = "PolymerNovelCompound"
    __table_args__ = (db.UniqueConstraint("name", "workbook"),)

    def __init__(self, **kwargs) -> None:
        self.time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    smiles = db.Column(db.Text, nullable=False)
    molec_weight = db.Column(db.Float(53))
    molec_formula = db.Column(db.Text, nullable=False)
    density = db.Column(db.Float(53))
    concentration = db.Column(db.Float(53))
    state = db.Column(db.Text)
    form = db.Column(db.Text)
    hphrase = db.Column(db.Text)
    workbook = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False, index=True
    )
    solvent = db.Column(db.ForeignKey("Solvent.name", ondelete="SET NULL"), index=True)
    time_of_creation = db.Column(db.DateTime)
