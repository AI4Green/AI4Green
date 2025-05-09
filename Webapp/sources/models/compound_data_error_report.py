from datetime import datetime

import pytz
from sources.extensions import db

from .base import Model


class CompoundDataErrorReport(Model):
    __tablename__ = "CompoundDataErrorReport"

    def __init__(self, **kwargs) -> None:
        self.time = datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        super().__init__(**kwargs)

    id = db.Column(db.Integer, primary_key=True)
    compound_name = db.Column(db.Text, nullable=False)
    compound = db.Column(
        db.ForeignKey("Compound.id", ondelete="CASCADE"), nullable=False, index=True
    )
    error_type = db.Column(db.Text, nullable=False)
    additional_info = db.Column(db.Text, nullable=False)
    time = db.Column(db.DateTime, nullable=False)
