from datetime import datetime
from enum import Enum

import pytz
from sources.extensions import db

from .base import Model

metadata = db.Model.metadata


class ReactionSet(Model):
    """
    A ReactionSet is a collection of related reactions run together as part of a high throughput screen or similar
    """

    __tablename__ = "ReactionSet"

    def __init__(self, **kwargs) -> None:
        self.time_of_creation = datetime.now(pytz.timezone("Europe/London")).replace(
            tzinfo=None
        )
        super().__init__(**kwargs)

    __table_args__ = (
        db.UniqueConstraint("workbook_id", "name"),
        db.UniqueConstraint("workbook_id", "set_id"),
    )

    id = db.Column(db.Integer, primary_key=True)
    set_id = db.Column(db.Text, nullable=False)

    reactions = db.relationship("Reaction", back_populates="reaction_set")

    name = db.Column(db.Text, nullable=False)
    description = db.Column(db.Text, nullable=False, default="")

    creator_id = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    creator = db.relationship("Person", foreign_keys=[creator_id])

    time_of_creation = db.Column(db.DateTime, nullable=False)
    time_of_update = db.Column(db.DateTime)
    date = db.Column(db.DateTime)

    reactor_dimensions = db.Column(db.JSON, nullable=False)

    workbook_id = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False
    )
    workbook = db.relationship("WorkBook", foreign_keys=[workbook_id])

    complete = db.Column(db.Text, nullable=False)
    status = db.Column(db.Text, nullable=False)

    # file_attachments = db.relationship("ReactionDataFile", backref="reaction_file")
    # addenda = db.relationship("ReactionNote", backref="reaction_addenda")

    data_export_request_id = db.Column(
        db.Integer, db.ForeignKey("DataExportRequest.id")
    )
