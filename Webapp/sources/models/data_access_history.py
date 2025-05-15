from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model


class DataAccessHistory(Model):
    __tablename__ = "DataAccessHistory"

    def __init__(self, **kwargs) -> None:
        self.time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)
    person_id = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    person = db.relationship("Person", backref="DataAccessHistory")
    time = db.Column(db.DateTime, nullable=False)
    workgroup_id = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )
    workbook_id = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=True, index=True
    )
    old_role = db.Column(db.Text, nullable=False)
    new_role = db.Column(db.Text, nullable=False)


# TODO: many functions for account deletion - only relevant if changes need to be visable on a workgroup/book not profile.
# delete_profile/routes

# TODO: make js to get blueprint
# TODO: make frontend / csv generation

# refactor changes>history
