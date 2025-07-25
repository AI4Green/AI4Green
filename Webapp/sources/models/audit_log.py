from sources.extensions import db

from .base import Model


class AuditLogEvent(Model):
    __tablename__ = "AuditLogEvent"
    __bind_key__ = "audit_log"

    id = db.Column(db.Integer, primary_key=True)
    event_time = db.Column(db.DateTime, nullable=False, index=True)
    event_type = db.Column(db.Text, nullable=False)
    message = db.Column(db.JSON, nullable=False)
