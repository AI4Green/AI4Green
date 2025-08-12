from sqlalchemy import Column, Integer, DateTime, Text, JSON, Index
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class AuditLogEvent(Base):
    __tablename__ = "AuditLogEvent"

    id = Column(Integer, primary_key=True)
    event_time = Column(DateTime, nullable=False)
    event_type = Column(Text, nullable=False)
    message = Column(JSON, nullable=False)

    __table_args__ = (Index("ix_audit_log_event_event_time", "event_time"),)
