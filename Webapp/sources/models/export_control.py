from .base import Model
from sources.extensions import db

class ExportControl(Model):
    """
    Database table to track usages of controlled chemical substances
    """
    __tablename__ = 'export_control'

    creator = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )

    workgroup = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=False
    )

    workbook = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False
    )

    reaction = db.Column(
        db.ForeignKey("Reaction.id", ondelete="CASCADE"), nullable=False
    )

    date_used = db.Column(db.DateTime)

    restricted_chemical_name = db.Column(db.Text, nullable=False)
    restricted_chemical_smiles = db.Column(db.Text, nullable=False)
    restricted_chemical_cas = db.Column(db.Text, nullable=False)

