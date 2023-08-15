from sources.extensions import db
from .base import Model


class WBStatusRequest(Model):
    __tablename__ = "WBStatusRequest"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    pi_sr = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    person = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    wb = db.Column(
        db.ForeignKey("WorkBook.id", ondelete="CASCADE"), nullable=False, index=True
    )
    current_role = db.Column(db.Text, nullable=False)
    new_role = db.Column(db.Text, nullable=False)
    time = db.Column(db.DateTime, nullable=False)
    notification = db.Column(
        db.ForeignKey("Notification.id", ondelete="CASCADE"), nullable=False, index=True
    )
    status = db.Column(db.Text, nullable=False)
    Person = db.relationship(
        "Person", backref="WBStatusRequestPerson", foreign_keys=[person]
    )
    Notification = db.relationship("Notification", backref="WBStatusRequest")
