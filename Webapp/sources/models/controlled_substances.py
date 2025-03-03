from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model


class ControlledSubstanceUsage(Model):
    """
    Database table to track usages of controlled chemical substances
    """

    __tablename__ = "controlled_substance_usage"

    def __init__(self, **kwargs) -> None:
        self.time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)

    creator = db.Column(db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False)
    creator_person = db.relationship("Person", foreign_keys=[creator])

    workgroups = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False
    )
    workgroup = db.relationship("WorkGroup", foreign_keys=[workgroups])

    workbooks = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False
    )
    workbook = db.relationship("WorkBook", foreign_keys=[workbooks])

    reactions = db.Column(
        db.ForeignKey("Reaction.id", ondelete="CASCADE"), nullable=False
    )
    reaction = db.relationship("Reaction", foreign_keys=[reactions])

    time_of_creation = db.Column(db.DateTime)

    last_edited = db.Column(db.DateTime)

    last_location = db.Column(db.Text, nullable=False)

    controlled_substance_name = db.Column(db.Text, nullable=True)

    controlled_substance_smiles = db.Column(db.Text, nullable=False)

    controlled_substance_cas = db.Column(db.Text, nullable=True)

    controlled_substance_inchi = db.Column(db.Text, nullable=False)
