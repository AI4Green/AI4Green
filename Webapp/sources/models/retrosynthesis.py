from datetime import datetime

import pytz
from flask_login import current_user
from sources.extensions import db

from .base import Model

metadata = db.Model.metadata


class Retrosynthesis(Model):
    """
    A Retrosynthesis is defined by;
    name: identifier
    date: date of creation of the reaction
    date_reaction: date when the reaction was performed
    precursor_reaction: links to the reaction in the database that preceded it
    successor_reaction: links to the reactions that follow on from this one
    green_metric: the score given for this reaction
    reactants: links to a set of reactants for this reaction
    products: links to a set of products for this reaction
    reagents: links to a set of reagents for this reaction
    solvent: links to a set of solvents for this reaction
    name_book_comp: this is the Unique PrimaryKey id for the reaction based on the
        name of the reaction and the workbook it is performed in
    """

    def __init__(self, **kwargs) -> None:
        self.time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        # self.creator = current_user.person
        super().__init__(**kwargs)

    __tablename__ = "Retrosynthesis"
    __table_args__ = (db.UniqueConstraint("workbook", "name"),)

    id = db.Column(db.Integer, primary_key=True)

    name = db.Column(db.Text, nullable=False)
    description = db.Column(db.Text, nullable=False, default="")
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
