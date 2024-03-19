from datetime import datetime, timedelta
from enum import Enum

import pytz
from sources.extensions import db

from .base import Model


class ApprovalStatus(Enum):
    PENDING = "pending"
    APPROVED = "approved"
    REJECTED = "rejected"
    EXPIRED = "expired"


class ExportFormat(Enum):
    RDF = "rdf"
    PDF = "pdf"
    ELN = "eln"
    SURF = "surf"
    CSV = "csv"


# Define an association table to represent the many-to-many relationship
data_export_request_approvers = db.Table(
    "data_export_request_approvers",
    db.Model.metadata,
    db.Column(
        "data_export_request_id", db.Integer, db.ForeignKey("DataExportRequest.id")
    ),
    db.Column("person_id", db.Integer, db.ForeignKey("Person.id")),
    db.Column("approved", db.Boolean, default=False),
)


class DataExportRequest(Model):
    """
    A Data export request is made for a workbook. It will export the contents in a user's chosen format
    """

    __tablename__ = "DataExportRequest"

    time_of_request = db.Column(
        db.DateTime,
        nullable=False,
        default=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None),
    )
    time_of_release = db.Column(
        db.DateTime,
        nullable=False,
        default=datetime.now(pytz.timezone("Europe/London")).replace(tzinfo=None)
        + timedelta(hours=48),
    )

    id = db.Column(db.Integer, primary_key=True)
    requestor = db.Column(
        db.ForeignKey("Person.id", ondelete="CASCADE"), nullable=False, index=True
    )
    requestor_person = db.relationship("Person", foreign_keys=[requestor])

    # Establish a many-to-many relationship with Person for required_approvers
    required_approvers = db.relationship(
        "Person",
        secondary=data_export_request_approvers,
        backref="data_export_requests",
    )

    status = db.Column(db.Enum(ApprovalStatus), default=ApprovalStatus.PENDING.value)

    workbook = db.Column(db.ForeignKey("WorkBook.id"), nullable=False)
    workbook = db.relationship("WorkBook", foreign_keys=[workbook])
    workgroup = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )
    WorkGroup = db.relationship(
        "WorkGroup", backref="WorkGroupRequest", cascade="all, delete"
    )
    institution = db.Column(
        db.ForeignKey("Institution.id", ondelete="CASCADE"), nullable=False, index=True
    )
    Institution = db.relationship("Institution", backref="WorkGroupRequest")
