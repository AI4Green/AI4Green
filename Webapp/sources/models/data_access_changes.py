from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model


class DataAccessChanges(Model):
    __tablename__ = "DataAccessChanges"

    def __init__(self, **kwargs) -> None:
        self.time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)
    person_id = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    person = db.relationship("Person", backref="DataAccessChanges")
    time = db.Column(db.DateTime, nullable=False)
    workgroup_id = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )
    old_role = db.Column(db.Text, nullable=False)
    new_role = db.Column(db.Text, nullable=False)


# TODO: write services to fill db on changes
# TODO: write services to access db
# TODO: check display string
# TODO: put function on account creation??? mightve already done this
