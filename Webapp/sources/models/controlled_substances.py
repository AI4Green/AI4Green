from .base import Model
from sources.extensions import db


class ControlledSubstanceUsage(Model):
    """
    Database table to track usages of controlled chemical substances
    """
    __tablename__ = 'controlled_substance_usage'

    id = db.Column(db.Integer, primary_key=True)

    creator = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False
    )

    workgroup = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False
    )

    workbook = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False
    )

    reaction = db.Column(
        db.ForeignKey("Reaction.id", ondelete="CASCADE"), nullable=False
    )

    date_created = db.Column(db.DateTime)

    last_edited = db.Column(db.DateTime)

    controlled_substance_name = db.Column(db.Text, nullable=False)

    controlled_substance_smiles = db.Column(db.Text, nullable=False)

    controlled_substance_cas = db.Column(db.Text, nullable=False)

