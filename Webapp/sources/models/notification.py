from sources.extensions import db
from .base import Model


class Notification(Model):
    __tablename__ = "Notification"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    person = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    type = db.Column(db.Text, nullable=False)
    info = db.Column(db.Text, nullable=False)
    time = db.Column(db.DateTime, nullable=False)
    status = db.Column(db.Text, nullable=False)
    wg = db.Column(db.Text, nullable=False, default="")
    wb = db.Column(db.Text, nullable=False, default="")
    wg_request = db.Column(db.Text, nullable=False, default="")
