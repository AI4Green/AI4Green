from sources.extensions import db
from .base import Model


class Compound(Model):
    __tablename__ = "Compound"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    cid = db.Column(db.Integer, unique=True)
    cas = db.Column(db.Text, nullable=False, unique=True)
    name = db.Column(db.Text, nullable=False)
    smiles = db.Column(db.Text)
    inchi = db.Column(db.Text)
    inchikey = db.Column(db.Text)
    molec_formula = db.Column(db.Text)
    density = db.Column(db.Float(53))
    concentration = db.Column(db.Float(53))
    boiling_point = db.Column(db.Float(53))
    melting_point = db.Column(db.Float(53))
    flash_point = db.Column(db.Float(53))
    autoignition_temp = db.Column(db.Float(53))
    molec_weight = db.Column(db.Float(53))
    state = db.Column(db.Text)
    form = db.Column(db.Text)
    hphrase = db.Column(db.Text)
    safety_score = db.Column(db.Float(53))
    health_score = db.Column(db.Float(53))
    enviro_score = db.Column(db.Float(53))
    econom_score = db.Column(db.Float(53))
    solvent = db.Column(db.ForeignKey("Solvent.name", ondelete="SET NULL"), index=True)
    error_report = db.relationship("CompoundDataErrorReport")
