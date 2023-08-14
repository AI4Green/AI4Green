from sources.extensions import db
from .base import Model


class WGStatusRequest(Model):
    __tablename__ = "WGStatusRequest"

    id = db.Column(
        db.Integer,
        primary_key=True
    )
    principal_investigator = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    person = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    wg = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )
    current_role = db.Column(db.Text, nullable=False)
    new_role = db.Column(db.Text, nullable=False)
    time = db.Column(db.DateTime, nullable=False)
    notification = db.Column(
        db.ForeignKey("Notification.id", ondelete="CASCADE"), nullable=False, index=True
    )
    status = db.Column(db.Text, nullable=False)
    Notification = db.relationship("Notification", backref="WGStatusRequest")
    workgroup = db.relationship("WorkGroup", backref="WGStatusRequest")
    pi = db.relationship(
        "Person", backref="WGStatusRequestPI", foreign_keys=[principal_investigator]
    )
    Person = db.relationship(
        "Person", backref="WGStatusRequestPerson", foreign_keys=[person]
    )
