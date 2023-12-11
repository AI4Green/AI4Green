from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model


class NovelCompound(Model):
    __tablename__ = "NovelCompound"
    __table_args__ = (
        db.UniqueConstraint("inchi", "workbook"),
        db.UniqueConstraint("name", "workbook"),
    )

    def __init__(self, **kwargs) -> None:
        self.time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.Text, nullable=False)
    cas = db.Column(db.Text)
    smiles = db.Column(db.Text, nullable=False)
    inchi = db.Column(db.Text)
    inchikey = db.Column(db.Text)
    molec_weight = db.Column(db.Float(53))
    molec_formula = db.Column(db.Text, nullable=False)
    density = db.Column(db.Float(53))
    concentration = db.Column(db.Float(53))
    boiling_point = db.Column(db.Float(53))
    melting_point = db.Column(db.Float(53))
    flash_point = db.Column(db.Float(53))
    autoignition_temp = db.Column(db.Float(53))
    state = db.Column(db.Text)
    form = db.Column(db.Text)
    hphrase = db.Column(db.Text)
    safety_score = db.Column(db.Float(53))
    health_score = db.Column(db.Float(53))
    enviro_score = db.Column(db.Float(53))
    econom_score = db.Column(db.Float(53))
    workbook = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False, index=True
    )
    solvent = db.Column(db.ForeignKey("Solvent.name", ondelete="SET NULL"), index=True)
    time_of_creation = db.Column(db.DateTime)
