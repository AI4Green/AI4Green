from sources.extensions import db
from .base import Model


class Solvent(Model):
    __tablename__ = "Solvent"

    name = db.Column(db.Text, primary_key=True)
    flag = db.Column(db.Integer)
    hazard = db.Column(db.Text)
    time_of_creation = db.Column(db.DateTime)

    compound = db.relationship("Compound")
    novel_compound = db.relationship("NovelCompound")
