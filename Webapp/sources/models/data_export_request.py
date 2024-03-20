from datetime import datetime, timedelta
from enum import Enum

import pytz
from sources.extensions import db

from .base import Model


class ApprovalStatus(Enum):
    PENDING = "PENDING"
    APPROVED = "APPROVED"
    REJECTED = "REJECTED"
    EXPIRED = "EXPIRED"


class ExportFormat(Enum):
    RDF = "RDF"
    PDF = "PDF"
    ELN = "ELN"
    SURF = "SURF"
    CSV = "CSV"
    SI = "SI"


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

# Association table for the many-to-many relationship between DataExportRequest and WorkBook
data_export_request_workbooks = db.Table(
    "data_export_request_workbook_association",
    db.Model.metadata,
    db.Column(
        "data_export_request_id", db.Integer, db.ForeignKey("DataExportRequest.id")
    ),
    db.Column("workbook_id", db.Integer, db.ForeignKey("WorkBook.id")),
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
        + timedelta(hours=168),
    )

    id = db.Column(db.Integer, primary_key=True)
    data_format = db.Column(db.Enum(ExportFormat))

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

    # workbooks_id = db.Column(db.Integer, db.ForeignKey("WorkBook.id"), nullable=False)
    # workbooks = db.relationship("WorkBook", foreign_keys=[workbooks_id])

    # reactions = db.relationship("Reaction", backref="DataExportRequest")
    reactions = db.relationship("Reaction", backref="data_export_request")

    # reactions = db.relationship("Reaction", secondary=)
    #
    # reactions = db.Column(
    #     db.ForeignKey("Reaction.id", ondelete="CASCADE"), nullable=False
    # )
    # Reactions = db.relationship("Reaction", foreign_keys=[reactions])

    # supports multiple workbooks if desired in future functionality
    workbooks = db.relationship(
        "WorkBook",
        secondary=data_export_request_workbooks,
        backref="data_export_requests",
    )

    workgroup = db.Column(
        db.ForeignKey("WorkGroup.id", ondelete="CASCADE"), nullable=False, index=True
    )
    WorkGroup = db.relationship(
        "WorkGroup", backref="data_export_request", cascade="all, delete"
    )
    institution = db.Column(
        db.ForeignKey("Institution.id", ondelete="CASCADE"), nullable=False, index=True
    )
    Institution = db.relationship("Institution", backref="data_export_request")
