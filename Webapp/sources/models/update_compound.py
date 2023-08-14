from typing import Dict
from sources.extensions import db
from .base import Model


class UpdateCompound(Model):
    """
    Update compound, specifically used in the 'Update' database.
    """
    __tablename__ = "UpdateCompound"
    __bind_key__ = 'update'

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


    def to_dict(self) -> Dict[str, str]:
        """
        Quick method to serialise the update to a dict.
        """
        return {
            'cid': self.cid,
            'cas': self.cas,
            'name': self.name,
            'smiles': self.smiles,
            'inchi': self.inchi,
            'inchikey': self.inchikey,
            'molec_formula': self.molec_formula,
            'density': self.density,
            'concentration': self.concentration,
            'boiling_point': self.boiling_point,
            'melting_point': self.melting_point,
            'flash_point': self.flash_point,
            'autoignition_temp': self.autoignition_temp,
            'molec_weight': self.molec_weight,
            'state': self.state,
            'form': self.form,
            'hphrase': self.hphrase,
            'safety_score': self.safety_score,
            'health_score': self.health_score,
            'enviro_score': self.enviro_score,
            'econom_score': self.econom_score,
        }